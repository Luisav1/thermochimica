#!/usr/bin/env python3
"""Generate committee-ready plots from the Fortran Hessian verification reports.

The Fortran executables remain the source of pass/fail decisions. This script
only runs their ``--report`` modes, preserves the complete step sweeps, and
plots scaled error against perturbation size with reference-order guides.
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit(
        "Matplotlib is required to generate Hessian verification figures.\n"
        "Install it in the Thermochimica container with:\n"
        "  apt-get update && apt-get install -y python3-matplotlib"
    ) from None

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
COLORS = ("#8a1538", "#006747", "#2f5597", "#c55a11", "#7030a0", "#444444")
MARKERS = ("o", "s", "^", "D", "v", "P")


def run_report(bin_dir: Path, executable: str) -> list[str]:
    bin_dir = bin_dir.resolve()
    result = subprocess.run(
        [str(bin_dir / executable), "--report"],
        cwd=bin_dir,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


def numbers(line: str) -> list[float]:
    return [float(value.replace("D", "E").replace("d", "e")) for value in re.findall(NUMBER, line)]


def numeric_rows(lines: list[str], header_index: int, minimum_columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in lines[header_index + 1 :]:
        values = numbers(line)
        if len(values) >= minimum_columns and line.lstrip()[:1].isdigit():
            rows.append(values)
        elif rows:
            break
    return rows


def fixed_table_rows(lines: list[str], header_index: int, columns: int) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in lines[header_index + 1 :]:
        fields = line.split()
        if len(fields) == columns:
            try:
                float(fields[0].replace("D", "E"))
            except ValueError:
                if rows:
                    break
                continue
            rows.append(fields)
        elif rows:
            break
    return rows


def find_after(lines: list[str], marker: str, header: str, start: int = 0) -> tuple[int, list[list[float]]]:
    marker_index = next(i for i in range(start, len(lines)) if marker in lines[i])
    header_index = next(i for i in range(marker_index + 1, len(lines)) if header in lines[i])
    return marker_index, numeric_rows(lines, header_index, 3)


def slope_guide(
    h: list[float], error: list[float], order: int, label: str
) -> tuple[str, list[float], list[float]] | None:
    valid = [(x, y) for x, y in zip(h, error) if x > 0.0 and y > 0.0 and math.isfinite(y)]
    if len(valid) < 2:
        return None
    anchor = min(2, len(valid) - 1)
    h0, e0 = valid[anchor]
    guide_h = [valid[0][0], valid[min(3, len(valid) - 1)][0]]
    # Deliberately offset the theoretical guide so it cannot be mistaken for a
    # duplicate measured series when the observed convergence is nearly exact.
    guide_e = [2.5 * e0 * (value / h0) ** order for value in guide_h]
    return label, guide_h, guide_e


def sweep_summary(
    label: str,
    h: list[float],
    error: list[float],
    order_min: float,
    order_max: float,
) -> tuple[str, str | None]:
    valid = [
        (step, value)
        for step, value in zip(h, error)
        if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
    ]
    if not valid:
        return f"{label}: no positive finite errors", None

    steps = [item[0] for item in valid]
    errors = [item[1] for item in valid]
    best_index = min(range(len(errors)), key=errors.__getitem__)
    in_range_orders: list[float] = []
    for index in range(best_index):
        if errors[index + 1] >= errors[index]:
            continue
        order = math.log(errors[index] / errors[index + 1]) / math.log(
            steps[index] / steps[index + 1]
        )
        if order_min <= order <= order_max:
            in_range_orders.append(order)

    if in_range_orders:
        order_text = (
            f"in-range observed $p$={min(in_range_orders):.2f}-{max(in_range_orders):.2f}"
        )
    else:
        order_text = "in-range observed $p$: N/A"

    roundoff_step = None
    for index in range(best_index + 1, len(errors)):
        if errors[index] > errors[index - 1]:
            roundoff_step = steps[index]
            break

    summary = (
        f"{label}: best={errors[best_index]:.2e} at $h$={steps[best_index]:.2e}; "
        f"{order_text}"
    )
    return summary, (
        f"{label} roundoff onset: $h$={roundoff_step:.2e}" if roundoff_step is not None else None
    )


def combined_sweep_summary(
    label: str,
    series: list[tuple[list[float], list[float]]],
    order_min: float,
    order_max: float,
) -> list[str]:
    best_errors: list[float] = []
    in_range_orders: list[float] = []
    roundoff_steps: list[float] = []
    for h, error in series:
        valid = [
            (step, value)
            for step, value in zip(h, error)
            if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
        ]
        if not valid:
            continue
        steps = [item[0] for item in valid]
        errors = [item[1] for item in valid]
        best_index = min(range(len(errors)), key=errors.__getitem__)
        best_errors.append(errors[best_index])
        for index in range(best_index):
            if errors[index + 1] >= errors[index]:
                continue
            order = math.log(errors[index] / errors[index + 1]) / math.log(
                steps[index] / steps[index + 1]
            )
            if order_min <= order <= order_max:
                in_range_orders.append(order)
        for index in range(best_index + 1, len(errors)):
            if errors[index] > errors[index - 1]:
                roundoff_steps.append(steps[index])
                break

    if not best_errors:
        return [f"{label}: no positive finite errors"]
    order_text = (
        f"in-range observed $p$={min(in_range_orders):.2f}-{max(in_range_orders):.2f}"
        if in_range_orders
        else "in-range observed $p$: N/A"
    )
    summary = [
        f"{label}: best-error range={min(best_errors):.2e}-{max(best_errors):.2e}; {order_text}"
    ]
    if roundoff_steps:
        if math.isclose(min(roundoff_steps), max(roundoff_steps), rel_tol=1e-12):
            summary.append(f"Roundoff onset across shown series: $h$={roundoff_steps[0]:.2e}")
        else:
            summary.append(
                f"Roundoff onset across shown series: $h$={min(roundoff_steps):.2e}-"
                f"{max(roundoff_steps):.2e}"
            )
    else:
        summary.append("No roundoff upturn in the retained sweep")
    return summary


def write_plot(
    output_stem: Path,
    title: str,
    series: list[tuple[str, list[float], list[float], bool]],
    subtitle: str | None = None,
    summary: list[str] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 6.7))
    plotted = 0
    for index, (label, xs, ys, dashed) in enumerate(series):
        usable = [
            (x, y)
            for x, y in zip(xs, ys)
            if x > 0.0 and y > 0.0 and math.isfinite(x) and math.isfinite(y)
        ]
        if len(usable) < 2:
            continue
        x_values, y_values = zip(*sorted(usable))
        color = COLORS[index % len(COLORS)]
        if dashed:
            ax.plot(x_values, y_values, linestyle="--", linewidth=1.8, color=color, label=label)
        else:
            ax.plot(
                x_values,
                y_values,
                marker=MARKERS[index % len(MARKERS)],
                markersize=5,
                linewidth=1.8,
                color=color,
                label=label,
            )
        plotted += 1

    if plotted == 0:
        raise ValueError(f"No positive finite data available for {title}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Perturbation size, $h$")
    ax.set_ylabel("Scale-normalized error")
    ax.set_title(title, fontweight="bold", pad=24 if subtitle else 12)
    if subtitle:
        ax.text(
            0.5,
            1.015,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=10,
            color="#444444",
        )
    ax.grid(which="major", color="#cfcfcf", linewidth=0.8)
    ax.grid(which="minor", color="#e8e8e8", linewidth=0.5, alpha=0.8)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    if summary:
        fig.text(
            0.5,
            0.018,
            "\n".join(summary),
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="#333333",
        )
    fig.subplots_adjust(left=0.13, right=0.97, top=0.90, bottom=0.31 if summary else 0.25)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rkmp(lines: list[str], output_dir: Path) -> None:
    marker = next(i for i, line in enumerate(lines) if "controlled energy direction" in line)
    header = next(i for i in range(marker, len(lines)) if "RKMP3 abs" in lines[i])
    rows = fixed_table_rows(lines, header, 7)
    h = [float(row[0]) for row in rows]
    err3 = [float(row[2]) for row in rows]
    err5 = [float(row[5]) for row in rows]
    plot_series = [
        ("RKMP controlled excess, 3-point", h, err3, False),
        ("RKMP controlled excess, 5-point", h, err5, False),
    ]
    summary3, roundoff3 = sweep_summary("3-point", h, err3, 1.7, 2.3)
    summary5, roundoff5 = sweep_summary("5-point", h, err5, 3.2, 4.8)
    summary = [summary3, summary5]
    if roundoff3 or roundoff5:
        summary.append("; ".join(item for item in (roundoff3, roundoff5) if item))
    for guide in (
        slope_guide(h, err3, 2, r"Reference slope $O(h^2)$"),
        slope_guide(h, err5, 4, r"Reference slope $O(h^4)$"),
    ):
        if guide:
            plot_series.append((*guide, True))
    write_plot(
        output_dir / "rkmp_energy",
        "RKMP controlled excess-energy verification",
        plot_series,
        summary=summary,
    )

    plot_series = []
    production_start = next(
        i
        for i, line in enumerate(lines)
        if "production partial-molar RKMP comparison" in line
    )
    cursor = production_start + 1
    while True:
        try:
            marker = next(
                i for i in range(cursor, len(lines)) if lines[i].startswith("direction: species")
            )
        except StopIteration:
            break
        header = next(i for i in range(marker + 1, len(lines)) if "norm abs" in lines[i])
        rows = numeric_rows(lines, header, 3)
        cursor = marker + 1
        label = lines[marker].removeprefix("direction: ").replace(" minus ", " - ")
        h = [row[0] for row in rows]
        error = [row[2] for row in rows]
        if any(value > 0.0 for value in error):
            plot_series.append((label, h, error, False))
    if plot_series:
        best_errors = [
            min(value for value in error if value > 0.0 and math.isfinite(value))
            for _, _, error, _ in plot_series
        ]
        # The native TestThermo30 polynomial begins at roundoff-level agreement,
        # so this plot demonstrates accuracy rather than a truncation-order slope.
        write_plot(
            output_dir / "rkmp_partial_molar",
            "RKMP production partial-molar verification",
            plot_series,
            "Unmodified TestThermo30 parameters; agreement is roundoff-limited from the coarsest step",
            [
                f"Best scaled error across nonzero directions: {min(best_errors):.2e}-"
                f"{max(best_errors):.2e}; measured order: N/A"
            ],
        )


def plot_cef(lines: list[str], output_dir: Path) -> None:
    marker = next(i for i, line in enumerate(lines) if line.strip() == "ALABANDITE controlled")
    header = next(i for i in range(marker, len(lines)) if "3pt abs" in lines[i])
    rows = fixed_table_rows(lines, header, 8)
    rows = [row for row in rows if int(row[0]) == 1]
    h = [float(row[1]) for row in rows]
    err3 = [float(row[3]) for row in rows]
    err5 = [float(row[6]) for row in rows]
    plot_series = [
        ("CEF ALABANDITE controlled direction 1, 3-point", h, err3, False),
        ("CEF ALABANDITE controlled direction 1, 5-point", h, err5, False),
    ]
    summary3, roundoff3 = sweep_summary("3-point", h, err3, 1.7, 2.3)
    summary5, roundoff5 = sweep_summary("5-point", h, err5, 3.2, 4.8)
    summary = [summary3, summary5]
    if roundoff3 or roundoff5:
        summary.append("; ".join(item for item in (roundoff3, roundoff5) if item))
    for guide in (
        slope_guide(h, err3, 2, r"Reference slope $O(h^2)$"),
        slope_guide(h, err5, 4, r"Reference slope $O(h^4)$"),
    ):
        if guide:
            plot_series.append((*guide, True))
    write_plot(
        output_dir / "cef_energy",
        "CEF controlled scalar-energy verification",
        plot_series,
        summary=summary,
    )


def plot_standalone_mqmqa(lines: list[str], output_dir: Path) -> None:
    plot_series = []
    guide_data: tuple[list[float], list[float]] | None = None
    summary_series: list[tuple[list[float], list[float]]] = []
    for case in ("G binary", "Q binary", "B family"):
        marker = next(i for i, line in enumerate(lines) if line.strip() == case)
        header = next(i for i in range(marker, len(lines)) if "3pt abs" in lines[i])
        rows = numeric_rows(lines, header, 7)
        h = [row[0] for row in rows]
        error = [row[5] for row in rows]
        plot_series.append((f"standalone MQMQA {case}, 5-point", h, error, False))
        summary_series.append((h, error))
        guide_data = guide_data or (h, error)
    if guide_data:
        guide = slope_guide(*guide_data, 4, r"Reference slope $O(h^4)$")
        if guide:
            plot_series.append((*guide, True))
    write_plot(
        output_dir / "mqmqa_standalone",
        "Standalone MQMQA G/Q/B verification",
        plot_series,
        summary=combined_sweep_summary("G/Q/B 5-point", summary_series, 3.2, 4.8),
    )


def plot_native_mqmqa(lines: list[str], output_dir: Path) -> None:
    header = next(i for i, line in enumerate(lines) if line.startswith("dir  h"))
    rows = []
    for line in lines[header + 1 :]:
        fields = line.split()
        if len(fields) != 8:
            continue
        try:
            rows.append([float(field) if field != "N/A" else math.nan for field in fields])
        except ValueError:
            continue
    plot_series = []
    guide_data: tuple[list[float], list[float]] | None = None
    summary_series: list[tuple[list[float], list[float]]] = []
    for direction in sorted({int(row[0]) for row in rows}):
        selected = [row for row in rows if int(row[0]) == direction]
        h = [row[1] for row in selected]
        error = [row[3] for row in selected]
        plot_series.append((f"native MQMQA direction {direction}", h, error, False))
        summary_series.append((h, error))
        guide_data = guide_data or (h, error)
    if guide_data:
        guide = slope_guide(*guide_data, 2, r"Reference slope $O(h^2)$")
        if guide:
            plot_series.append((*guide, True))
    write_plot(
        output_dir / "mqmqa_native",
        "Native MQMQA G-family verification",
        plot_series,
        summary=combined_sweep_summary("Five native directions", summary_series, 1.7, 2.3),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=Path, default=Path("bin"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/hessian_verification"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reports = {
        "rkmp": run_report(args.bin_dir, "TestRKMPHessianVerification"),
        "cef": run_report(args.bin_dir, "TestCEFHessianVerification"),
        "standalone": run_report(args.bin_dir, "TestMQMQAHessianVerification"),
        "native": run_report(args.bin_dir, "TestMQMQANativeHessianVerification"),
    }
    for name, lines in reports.items():
        (args.output_dir / f"{name}_report.txt").write_text("\n".join(lines) + "\n")

    plot_rkmp(reports["rkmp"], args.output_dir)
    plot_cef(reports["cef"], args.output_dir)
    plot_standalone_mqmqa(reports["standalone"], args.output_dir)
    plot_native_mqmqa(reports["native"], args.output_dir)
    print(f"Wrote reports and plots to {args.output_dir}")


if __name__ == "__main__":
    main()
