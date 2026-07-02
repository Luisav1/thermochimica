#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

docker start thermochimica >/dev/null

docker exec -it thermochimica bash -lc "
  cd /work &&
  chmod +x scripts/watch-build.sh &&
  ./scripts/watch-build.sh
"