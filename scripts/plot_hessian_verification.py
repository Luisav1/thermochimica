#!/usr/bin/env python3
"""Generate committee-ready plots from the Fortran Hessian verification reports.

The Fortran executables remain the source of pass/fail decisions. This script
only runs their ``--report`` modes, preserves the complete step sweeps, and
plots scaled error against perturbation size with reference-order guides.
Displayed ``in-range observed p`` values are descriptive individual slopes;
unlike the Fortran gate, this plotting summary does not require two consecutive
acceptable orders with three consecutively decreasing errors.
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
    h: list[float],
    error: list[float],
    order: int,
    label: str,
    order_min: float | None = None,
    order_max: float | None = None,
) -> tuple[str, list[float], list[float]] | None:
    valid = [(x, y) for x, y in zip(h, error) if x > 0.0 and y > 0.0 and math.isfinite(y)]
    if len(valid) < 2:
        return None

    best_index = min(range(len(valid)), key=lambda index: valid[index][1])
    guide_start = 0
    guide_end = best_index if best_index > 0 else len(valid) - 1

    # When an expected-order window is supplied, use the longest consecutive
    # run of in-range, decreasing errors. This keeps the guide on the measured
    # truncation region instead of extending it into a coarse nonlinear point.
    if (
        order_min is not None
        and order_max is not None
        and best_index >= 1
        and best_index < len(valid) - 1
    ):
        in_range: list[bool] = []
        for index in range(best_index):
            h_coarse, error_coarse = valid[index]
            h_fine, error_fine = valid[index + 1]
            if error_fine >= error_coarse:
                in_range.append(False)
                continue
            observed = math.log(error_coarse / error_fine) / math.log(h_coarse / h_fine)
            in_range.append(order_min <= observed <= order_max)

        longest_start = 0
        longest_length = 0
        current_start = 0
        current_length = 0
        for index, accepted in enumerate(in_range):
            if accepted:
                if current_length == 0:
                    current_start = index
                current_length += 1
                if current_length > longest_length:
                    longest_start = current_start
                    longest_length = current_length
            else:
                current_length = 0
        if longest_length > 0:
            guide_start = longest_start
            guide_end = longest_start + longest_length

    guide_region = valid[guide_start : guide_end + 1]
    if len(guide_region) < 2:
        return None
    anchor = min(2, len(guide_region) - 1)
    h0, e0 = guide_region[anchor]
    guide_h = [value[0] for value in guide_region]
    # Deliberately offset the theoretical guide so it cannot be mistaken for a
    # duplicate measured series when the observed convergence is nearly exact.
    guide_e = [2.5 * e0 * (value / h0) ** order for value in guide_h]
    return label, guide_h, guide_e


def roundoff_onset(h: list[float], error: list[float]) -> float | None:
    valid = [
        (step, value)
        for step, value in zip(h, error)
        if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
    ]
    if len(valid) < 2:
        return None
    best_index = min(range(len(valid)), key=lambda index: valid[index][1])
    for index in range(best_index + 1, len(valid)):
        if valid[index][1] > valid[index - 1][1]:
            return valid[index][0]
    return None


def in_range_order_bounds(
    h: list[float], error: list[float], order_min: float, order_max: float
) -> tuple[float, float] | None:
    valid = [
        (step, value)
        for step, value in zip(h, error)
        if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
    ]
    if len(valid) < 2:
        return None
    best_index = min(range(len(valid)), key=lambda index: valid[index][1])
    orders: list[float] = []
    for index in range(best_index):
        h_coarse, error_coarse = valid[index]
        h_fine, error_fine = valid[index + 1]
        if error_fine >= error_coarse:
            continue
        observed = math.log(error_coarse / error_fine) / math.log(h_coarse / h_fine)
        if order_min <= observed <= order_max:
            orders.append(observed)
    return (min(orders), max(orders)) if orders else None


def set_measured_limits(
    ax: plt.Axes, series: list[tuple[list[float], list[float]]]
) -> None:
    measured = [
        (step, value)
        for h, error in series
        for step, value in zip(h, error)
        if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
    ]
    if not measured:
        return
    log_x = [math.log10(step) for step, _ in measured]
    log_y = [math.log10(value) for _, value in measured]
    x_pad = max(0.08, 0.03 * (max(log_x) - min(log_x)))
    y_pad = max(0.15, 0.04 * (max(log_y) - min(log_y)))
    ax.set_xlim(10 ** (min(log_x) - x_pad), 10 ** (max(log_x) + x_pad))
    ax.set_ylim(10 ** (min(log_y) - y_pad), 10 ** (max(log_y) + y_pad))


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
    # This is intentionally a presentation summary, not a second implementation
    # of AssessFDSweep. Report every individually in-range pre-minimum slope and
    # leave the stricter consecutive-order pass/fail decision to Fortran.
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
    caption: str | None = None,
    roundoff_start: float | None = None,
    roundoff_label: str = "Roundoff-dominated region",
) -> None:
    fig, ax = plt.subplots(figsize=(9.6, 6.8))
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
    set_measured_limits(
        ax,
        [(xs, ys) for _, xs, ys, dashed in series if not dashed],
    )
    ax.set_xlabel("Perturbation size, $h$ (smaller $h$ to the left)")
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

    if roundoff_start is not None:
        left_limit = ax.get_xlim()[0]
        right_limit = ax.get_xlim()[1]
        shade_right = min(max(roundoff_start, left_limit), right_limit)
        ax.axvspan(
            left_limit,
            shade_right,
            color="#777777",
            alpha=0.08,
            zorder=0,
            label=roundoff_label,
        )

    ax.legend(
        loc="upper left",
        ncol=1,
        frameon=False,
        fontsize=9,
    )
    if summary:
        ax.text(
            0.97,
            0.03,
            "\n".join(summary),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.88},
        )
    if caption:
        fig.text(
            0.5,
            0.025,
            caption,
            ha="center",
            va="bottom",
            fontsize=9,
            color="#333333",
            wrap=True,
        )
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.15 if caption else 0.11)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_rkmp_order_plot(
    output_stem: Path,
    h_one_sided: list[float],
    err_forward: list[float],
    err_backward: list[float],
    h_centered: list[float],
    err3: list[float],
    err5: list[float],
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.5))
    panels = (
        (
            axes[0],
            "A. One-sided second differences",
            [
                ("Forward", h_one_sided, err_forward, COLORS[0], MARKERS[0]),
                ("Backward", h_one_sided, err_backward, COLORS[1], MARKERS[1]),
            ],
            [
                slope_guide(
                    h_one_sided,
                    err_forward,
                    1,
                    r"Reference $O(h)$",
                    0.7,
                    1.3,
                )
            ],
            [
                ("forward", in_range_order_bounds(h_one_sided, err_forward, 0.7, 1.3)),
                ("backward", in_range_order_bounds(h_one_sided, err_backward, 0.7, 1.3)),
            ],
            [
                ("forward", roundoff_onset(h_one_sided, err_forward)),
                ("backward", roundoff_onset(h_one_sided, err_backward)),
            ],
        ),
        (
            axes[1],
            "B. Centered second differences",
            [
                ("Three-point", h_centered, err3, COLORS[2], MARKERS[2]),
                ("Five-point", h_centered, err5, COLORS[3], MARKERS[3]),
            ],
            [
                slope_guide(
                    h_centered,
                    err3,
                    2,
                    r"Reference $O(h^2)$",
                    1.7,
                    2.3,
                ),
                slope_guide(
                    h_centered,
                    err5,
                    4,
                    r"Reference $O(h^4)$",
                    3.2,
                    4.8,
                ),
            ],
            [
                ("three-point", in_range_order_bounds(h_centered, err3, 1.7, 2.3)),
                ("five-point", in_range_order_bounds(h_centered, err5, 3.2, 4.8)),
            ],
            [
                ("three-point", roundoff_onset(h_centered, err3)),
                ("five-point", roundoff_onset(h_centered, err5)),
            ],
        ),
    )

    for ax, title, measured, guides, order_bounds, onsets in panels:
        for label, h_values, error_values, color, marker in measured:
            usable = sorted(
                (step, value)
                for step, value in zip(h_values, error_values)
                if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
            )
            ax.plot(
                [item[0] for item in usable],
                [item[1] for item in usable],
                color=color,
                marker=marker,
                markersize=5,
                linewidth=1.8,
                label=label,
            )
        for guide_index, guide in enumerate(item for item in guides if item is not None):
            label, guide_h, guide_error = guide
            ax.plot(
                sorted(guide_h),
                [value for _, value in sorted(zip(guide_h, guide_error))],
                linestyle="--",
                linewidth=1.8,
                color=COLORS[4 + guide_index],
                label=label,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        set_measured_limits(ax, [(item[1], item[2]) for item in measured])
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Perturbation size, $h$\n(smaller $h$ to the left)")
        ax.set_ylabel("Scale-normalized error")
        ax.grid(which="major", color="#cfcfcf", linewidth=0.8)
        ax.grid(which="minor", color="#e8e8e8", linewidth=0.5, alpha=0.8)

        finite_onsets = [value for _, value in onsets if value is not None]
        if finite_onsets:
            left_limit = ax.get_xlim()[0]
            ax.axvspan(
                left_limit,
                max(finite_onsets),
                color="#777777",
                alpha=0.08,
                zorder=0,
            )

        order_lines = [
            f"{name} p={bounds[0]:.2f}-{bounds[1]:.2f}"
            for name, bounds in order_bounds
            if bounds is not None
        ]
        onset_lines = [
            f"{name} roundoff h={value:.2e}"
            for name, value in onsets
            if value is not None
        ]
        if not onset_lines:
            onset_lines = ["no roundoff upturn in retained sweep"]
        ax.text(
            0.97,
            0.03,
            "\n".join(order_lines + onset_lines),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.88},
        )
        ax.legend(loc="upper left", frameon=False, fontsize=9)

    fig.suptitle("RKMP controlled excess-energy verification", fontsize=17, fontweight="bold")
    fig.text(
        0.5,
        0.90,
        "Dashed lines show theoretical truncation rates over the measured convergence region.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.text(
        0.5,
        0.035,
        "Forward/backward absolute errors may overlap because their leading signed errors have "
        "opposite signs but similar magnitudes. Gray shading marks where roundoff affects at least one curve.",
        ha="center",
        fontsize=9,
        color="#333333",
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.83, bottom=0.16, wspace=0.25)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rkmp(lines: list[str], output_dir: Path) -> None:
    marker = next(i for i, line in enumerate(lines) if "controlled energy direction" in line)
    one_sided_header = next(i for i in range(marker, len(lines)) if "forward abs" in lines[i])
    one_sided_rows = fixed_table_rows(lines, one_sided_header, 7)
    h_one_sided = [float(row[0]) for row in one_sided_rows]
    err_forward = [float(row[2]) for row in one_sided_rows]
    err_backward = [float(row[5]) for row in one_sided_rows]
    header = next(i for i in range(marker, len(lines)) if "RKMP3 abs" in lines[i])
    rows = fixed_table_rows(lines, header, 7)
    h = [float(row[0]) for row in rows]
    err3 = [float(row[2]) for row in rows]
    err5 = [float(row[5]) for row in rows]
    write_rkmp_order_plot(
        output_dir / "rkmp_energy",
        h_one_sided,
        err_forward,
        err_backward,
        h,
        err3,
        err5,
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
            "Central finite differences of established production partial molars at the "
            "converged TestThermo30 state",
            [
                f"best scaled error={min(best_errors):.2e}-{max(best_errors):.2e}",
                "measured order: N/A",
                "roundoff-limited from coarsest retained step",
            ],
            "The unmodified low-order RKMP polynomial agrees at nearly machine precision "
            "before refinement. This is a native Thermochimica accuracy cross-check, not "
            "a truncation-order demonstration.",
            roundoff_start=max(
                step
                for _, h_values, _, _ in plot_series
                for step in h_values
                if step > 0.0 and math.isfinite(step)
            ),
            roundoff_label="Roundoff-limited sweep",
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
        ("Three-point", h, err3, False),
        ("Five-point", h, err5, False),
    ]
    onset3 = roundoff_onset(h, err3)
    onset5 = roundoff_onset(h, err5)
    onset_values = [value for value in (onset3, onset5) if value is not None]
    order3 = in_range_order_bounds(h, err3, 1.7, 2.3)
    order5 = in_range_order_bounds(h, err5, 3.2, 4.8)
    best3 = min(value for value in err3 if value > 0.0 and math.isfinite(value))
    best5 = min(value for value in err5 if value > 0.0 and math.isfinite(value))
    summary = [
        f"3-point: p={order3[0]:.2f}-{order3[1]:.2f}; best={best3:.2e}",
        f"5-point: p={order5[0]:.2f}-{order5[1]:.2f}; best={best5:.2e}",
        f"roundoff onset h={onset3:.2e} / {onset5:.2e}",
    ]
    for guide in (
        slope_guide(h, err3, 2, r"Reference $O(h^2)$", 1.7, 2.3),
        slope_guide(h, err5, 4, r"Reference $O(h^4)$", 3.2, 4.8),
    ):
        if guide:
            plot_series.append((*guide, True))
    write_plot(
        output_dir / "cef_energy",
        "CEF controlled scalar-energy verification",
        plot_series,
        "Three- and five-point directional curvatures for a positive-interior "
        "plain-SUBL ALABANDITE state",
        summary=summary,
        caption="The measured second- and fourth-order regions verify the standalone CEF "
        "Hessian against scalar-energy finite differences. Gray shading marks the "
        "identified small-h roundoff region.",
        roundoff_start=max(onset_values) if onset_values else None,
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
        plot_series.append((case, h, error, False))
        summary_series.append((h, error))
        guide_data = guide_data or (h, error)
    if guide_data:
        guide = slope_guide(*guide_data, 4, r"Reference $O(h^4)$", 3.2, 4.8)
        if guide:
            plot_series.append((*guide, True))
    onset_values = [
        onset
        for h_values, error_values in summary_series
        if (onset := roundoff_onset(h_values, error_values)) is not None
    ]
    order_values = [
        bounds
        for h_values, error_values in summary_series
        if (bounds := in_range_order_bounds(h_values, error_values, 3.2, 4.8)) is not None
    ]
    best_values = [
        min(value for value in error_values if value > 0.0 and math.isfinite(value))
        for _, error_values in summary_series
    ]
    if math.isclose(min(onset_values), max(onset_values), rel_tol=1e-12):
        onset_summary = f"roundoff onset h={onset_values[0]:.2e}"
    else:
        onset_summary = (
            f"roundoff onset h={min(onset_values):.2e}-{max(onset_values):.2e}"
        )
    write_plot(
        output_dir / "mqmqa_standalone",
        "Standalone MQMQA scalar-energy verification",
        plot_series,
        "Controlled nonmagnetic SUBG G-, Q-, and B-family cases using five-point "
        "directional differences",
        summary=[
            f"5-point p={min(item[0] for item in order_values):.2f}-"
            f"{max(item[1] for item in order_values):.2f}",
            f"best scaled error={min(best_values):.2e}-{max(best_values):.2e}",
            onset_summary,
        ],
        caption="The controlled cases independently check the traced standalone G, Q, and "
        "B scalar forms and their analytic Hessian propagation. They are not "
        "database-backed native Thermochimica comparisons.",
        roundoff_start=max(onset_values) if onset_values else None,
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
        plot_series.append((f"Direction {direction}", h, error, False))
        summary_series.append((h, error))
        guide_data = guide_data or (h, error)
    if guide_data:
        guide = slope_guide(*guide_data, 2, r"Reference $O(h^2)$", 1.7, 2.3)
        if guide:
            plot_series.append((*guide, True))
    order_values = [
        bounds
        for h_values, error_values in summary_series
        if (bounds := in_range_order_bounds(h_values, error_values, 1.7, 2.3)) is not None
    ]
    best_values = [
        min(value for value in error_values if value > 0.0 and math.isfinite(value))
        for _, error_values in summary_series
    ]
    write_plot(
        output_dir / "mqmqa_native",
        "Native MQMQA G-family verification",
        plot_series,
        "Central finite differences of production partial molars along five "
        "total-preserving directions",
        summary=[
            f"five directions: p={min(item[0] for item in order_values):.2f}-"
            f"{max(item[1] for item in order_values):.2f}",
            f"best scaled error={min(best_values):.2e}-{max(best_values):.2e}",
            "no roundoff upturn in retained sweep",
        ],
        caption="At a converged CuFeC-Kang.dat plain-SUBG state, analytic Hessian-vector "
        "products agree with production Thermochimica partial-molar differences. "
        "Scope: nonmagnetic reference/configurational/G-family behavior.",
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
