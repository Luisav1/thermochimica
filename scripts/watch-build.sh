#!/usr/bin/env bash

set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR" || exit 1

usage() {
  cat <<'EOF'
Usage: scripts/watch-build.sh [make target...]

Examples:
  scripts/watch-build.sh
  scripts/watch-build.sh obj/InitGEMSolver.o obj/MapRKMPHessianToGEMVariables.o
  scripts/watch-build.sh -j

Notes:
  - With no arguments, the script runs: make -j
  - If `watchexec`, `entr`, `inotifywait`, or `fswatch` is installed later,
    this script will use it automatically.
  - Otherwise it falls back to a portable polling loop.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$#" -eq 0 ]]; then
  MAKE_ARGS=(-j)
else
  MAKE_ARGS=("$@")
fi

# By default, run tests after a successful build. Disable with `-n` or
# `--no-run-tests`, or control with the `RUN_TESTS` environment variable.
RUN_TESTS=${RUN_TESTS:-1}
filtered=()
for _arg in "${MAKE_ARGS[@]}"; do
  case "${_arg}" in
    -T|--run-tests)
      RUN_TESTS=1
      ;;
    -n|--no-run-tests)
      RUN_TESTS=0
      ;;
    *)
      filtered+=("${_arg}")
      ;;
  esac
done
MAKE_ARGS=("${filtered[@]}")

# Build command string for use in watcher helpers that run a shell command
MAKE_CMD="make"
for _arg in "${MAKE_ARGS[@]}"; do
  MAKE_CMD+=" $(printf '%q' "${_arg}")"
done

run_build() {
  printf '\n[%s] Running: make' "$(date '+%Y-%m-%d %H:%M:%S')"
  for arg in "${MAKE_ARGS[@]}"; do
    printf ' %q' "$arg"
  done
  printf '\n\n'
  if make "${MAKE_ARGS[@]}"; then
    if [[ "${RUN_TESTS}" -eq 1 ]]; then
      run_tests
    fi
    return 0
  else
    return 1
  fi
}

run_tests() {
  printf '\n[%s] Running tests\n' "$(date '+%Y-%m-%d %H:%M:%S')"
  if [[ -x ./run_tests ]]; then
    ./run_tests
    return $?
  fi
  make test
  return $?
}

source_files() {
  find src test \
    \( -name '*.f90' -o -name '*.F90' -o -name '*.inc' -o -name '*.h' \) \
    -type f | sort
}

poll_signature() {
  source_files | xargs -r stat -c '%Y %n' 2>/dev/null
}

if command -v watchexec >/dev/null 2>&1; then
  if [[ "${RUN_TESTS}" -eq 1 ]]; then
    exec watchexec -e f90,F90,inc,h -- sh -c "${MAKE_CMD} && (./run_tests || make test)"
  else
    exec watchexec -e f90,F90,inc,h -- sh -c "${MAKE_CMD}"
  fi
fi

if command -v entr >/dev/null 2>&1; then
  if [[ "${RUN_TESTS}" -eq 1 ]]; then
    source_files | entr -r sh -c "${MAKE_CMD} && (./run_tests || make test)"
  else
    source_files | entr -r make "${MAKE_ARGS[@]}"
  fi
  exit $?
fi

if command -v inotifywait >/dev/null 2>&1; then
  run_build
  while inotifywait -qq -r -e close_write,create,delete,move src test; do
    run_build
  done
  exit $?
fi

if command -v fswatch >/dev/null 2>&1; then
  run_build
  fswatch -0 src test | while IFS= read -r -d '' _event; do
    run_build
  done
  exit $?
fi

echo "No file-watcher tool found; using polling mode (checks every 1s)."
echo "Press Ctrl+C to stop."

last_sig="$(poll_signature)"
run_build

while true; do
  sleep 1
  new_sig="$(poll_signature)"
  if [[ "$new_sig" != "$last_sig" ]]; then
    last_sig="$new_sig"
    run_build
  fi
done
