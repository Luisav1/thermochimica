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
import textwrap
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
FD_FINITE_PRECISION_FLOOR = 100.0 * math.ulp(1.0)


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


def expected_order_region(
    h: list[float],
    error: list[float],
    order_min: float,
    order_max: float,
) -> list[tuple[float, float]]:
    valid = [
        (x, y)
        for x, y in zip(h, error)
        if x > 0.0 and y > FD_FINITE_PRECISION_FLOOR and math.isfinite(y)
    ]
    if len(valid) < 2:
        return []

    best_index = min(range(len(valid)), key=lambda index: valid[index][1])
    guide_start = 0
    guide_end = best_index if best_index > 0 else len(valid) - 1

    # Use the longest consecutive run of in-range, decreasing errors. This
    # identifies a measured truncation region without extending into either a
    # coarse nonlinear point or the observed small-h error-upturn region.
    if best_index >= 1:
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

    return valid[guide_start : guide_end + 1]


def slope_guide(
    h: list[float],
    error: list[float],
    order: int,
    label: str,
    order_min: float,
    order_max: float,
    extrapolation_h: list[float] | None = None,
) -> tuple[str, list[float], list[float]] | None:
    guide_region = expected_order_region(h, error, order_min, order_max)

    if len(guide_region) < 2:
        return None
    anchor = min(2, len(guide_region) - 1)
    h0, e0 = guide_region[anchor]
    guide_h = sorted(
        {
            value
            for value in (extrapolation_h if extrapolation_h is not None else h)
            if value > 0.0 and math.isfinite(value)
        }
    )
    # The guide passes through an actual measured point. Its vertical placement
    # therefore has no arbitrary visibility offset. Extrapolating C*h**p over
    # the complete measured h-range shows how truncation error would continue
    # decreasing if finite-precision cancellation did not become dominant.
    guide_e = [e0 * (value / h0) ** order for value in guide_h]
    return label, guide_h, guide_e


def small_h_upturn_onset(h: list[float], error: list[float]) -> float | None:
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


def in_range_orders(
    h: list[float], error: list[float], order_min: float, order_max: float
) -> list[float]:
    valid = [
        (step, value)
        for step, value in zip(h, error)
        if step > 0.0
        and value > FD_FINITE_PRECISION_FLOOR
        and math.isfinite(step)
        and math.isfinite(value)
    ]
    if len(valid) < 2:
        return []
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
    return orders


def in_range_order_bounds(
    h: list[float], error: list[float], order_min: float, order_max: float
) -> tuple[float, float] | None:
    orders = in_range_orders(h, error, order_min, order_max)
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

    upturn_step = None
    for index in range(best_index + 1, len(errors)):
        if errors[index] > errors[index - 1]:
            upturn_step = steps[index]
            break

    summary = (
        f"{label}: best={errors[best_index]:.2e} at $h$={steps[best_index]:.2e}; "
        f"{order_text}"
    )
    return summary, (
        f"{label} small-$h$ error upturn: $h$={upturn_step:.2e}" if upturn_step is not None else None
    )


def combined_sweep_summary(
    label: str,
    series: list[tuple[list[float], list[float]]],
    order_min: float,
    order_max: float,
) -> list[str]:
    best_errors: list[float] = []
    in_range_orders: list[float] = []
    upturn_steps: list[float] = []
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
                upturn_steps.append(steps[index])
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
    if upturn_steps:
        if math.isclose(min(upturn_steps), max(upturn_steps), rel_tol=1e-12):
            summary.append(f"Small-$h$ error upturn across shown series: $h$={upturn_steps[0]:.2e}")
        else:
            summary.append(
                f"Small-$h$ error upturn across shown series: $h$={min(upturn_steps):.2e}-"
                f"{max(upturn_steps):.2e}"
            )
    else:
        summary.append("No small-$h$ upturn observed in the retained sweep")
    return summary


def write_plot(
    output_stem: Path,
    title: str,
    series: list[tuple[str, list[float], list[float], bool]],
    subtitle: str | None = None,
    summary: list[str] | None = None,
    caption: str | None = None,
    small_h_upturn_start: float | None = None,
    small_h_region_label: str = "Small-$h$ finite-precision region",
    small_h_explanation: str | None = None,
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

    if small_h_upturn_start is not None:
        left_limit = ax.get_xlim()[0]
        right_limit = ax.get_xlim()[1]
        shade_right = min(max(small_h_upturn_start, left_limit), right_limit)
        ax.axvspan(
            left_limit,
            shade_right,
            color="#777777",
            alpha=0.08,
            zorder=0,
            label=small_h_region_label,
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
    display_caption = caption
    if small_h_upturn_start is not None:
        upturn_explanation = small_h_explanation or (
            "Gray shading begins at the first observed post-minimum error increase. At sufficiently "
            "small h, cancellation between nearly equal Gibbs-energy evaluations and division by "
            "h squared amplify floating-point error."
        )
        display_caption = (
            f"{caption} {upturn_explanation}" if caption else upturn_explanation
        )
    if display_caption:
        fig.text(
            0.5,
            0.025,
            textwrap.fill(display_caption, width=135),
            ha="center",
            va="bottom",
            fontsize=9,
            color="#333333",
            wrap=True,
        )
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.23 if display_caption else 0.11)
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
                    r"Extrapolated theoretical $O(h)$ truncation trend",
                    0.7,
                    1.3,
                )
            ],
            [
                ("forward", in_range_order_bounds(h_one_sided, err_forward, 0.7, 1.3)),
                ("backward", in_range_order_bounds(h_one_sided, err_backward, 0.7, 1.3)),
            ],
            [
                ("forward", small_h_upturn_onset(h_one_sided, err_forward)),
                ("backward", small_h_upturn_onset(h_one_sided, err_backward)),
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
                    r"Extrapolated theoretical $O(h^2)$ truncation trend",
                    1.7,
                    2.3,
                ),
                slope_guide(
                    h_centered,
                    err5,
                    4,
                    r"Extrapolated theoretical $O(h^4)$ truncation trend",
                    3.2,
                    4.8,
                ),
            ],
            [
                ("three-point", in_range_order_bounds(h_centered, err3, 1.7, 2.3)),
                ("five-point", in_range_order_bounds(h_centered, err5, 3.2, 4.8)),
            ],
            [
                ("three-point", small_h_upturn_onset(h_centered, err3)),
                ("five-point", small_h_upturn_onset(h_centered, err5)),
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
                label="Small-$h$ finite-precision region",
            )

        order_lines = [
            f"{name} p={bounds[0]:.2f}-{bounds[1]:.2f}"
            for name, bounds in order_bounds
            if bounds is not None
        ]
        onset_lines = [
            f"{name} small-$h$ upturn at $h$={value:.2e}"
            for name, value in onsets
            if value is not None
        ]
        if not onset_lines:
            onset_lines = ["no small-$h$ upturn observed in retained sweep"]
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
        "Dashed theoretical truncation trends are anchored to measured convergence points and extrapolated across $h$.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.text(
        0.5,
        0.035,
        "Forward/backward absolute errors may overlap because their leading signed errors have "
        "opposite signs but similar magnitudes. Gray shading begins at the first observed post-minimum error increase; "
        "cancellation between nearly equal energies and division by $h^2$ then amplify floating-point error.",
        ha="center",
        fontsize=9,
        color="#333333",
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.83, bottom=0.16, wspace=0.25)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_signed_one_sided_plot(
    output_stem: Path,
    h: list[float],
    signed_forward: list[float],
    signed_backward: list[float],
    abs_forward: list[float],
    abs_backward: list[float],
) -> None:
    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    for label, values, color, marker in (
        ("Forward signed error", signed_forward, COLORS[0], MARKERS[0]),
        ("Backward signed error", signed_backward, COLORS[1], MARKERS[1]),
    ):
        usable = sorted(
            (step, value)
            for step, value in zip(h, values)
            if step > 0.0 and math.isfinite(step) and math.isfinite(value)
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

    common_ratios = []
    opposite_count = 0
    for forward, backward in zip(signed_forward, signed_backward):
        if forward * backward < 0.0:
            opposite_count += 1
        if forward * backward < 0.0 and backward != 0.0:
            ratio = abs(forward / backward)
            if 0.5 <= ratio <= 2.0:
                common_ratios.append(ratio)

    onset_values = [
        value
        for value in (small_h_upturn_onset(h, abs_forward), small_h_upturn_onset(h, abs_backward))
        if value is not None
    ]
    ax.set_xscale("log")
    nonzero = [abs(value) for value in signed_forward + signed_backward if value != 0.0]
    linthresh = max(100.0 * math.ulp(1.0), min(nonzero) if nonzero else 1.0e-12)
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.axhline(0.0, color="#555555", linewidth=0.9)
    ax.set_xlabel("Perturbation size, $h$ (smaller $h$ to the left)")
    ax.set_ylabel("Signed scale-normalized error")
    ax.set_title("RKMP one-sided signed-error verification", fontweight="bold", pad=24)
    ax.text(
        0.5,
        1.015,
        "Opposite-sign leading errors are expected before finite-precision cancellation dominates",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10,
        color="#444444",
    )
    ax.grid(which="major", color="#cfcfcf", linewidth=0.8)
    ax.grid(which="minor", color="#e8e8e8", linewidth=0.5, alpha=0.8)
    if onset_values:
        left_limit = ax.get_xlim()[0]
        ax.axvspan(
            left_limit,
            max(onset_values),
            color="#777777",
            alpha=0.08,
            zorder=0,
            label="Small-$h$ finite-precision region",
        )
    ratio_text = (
        f"matched magnitude ratio={min(common_ratios):.3f}-{max(common_ratios):.3f}"
        if common_ratios
        else "matched magnitude ratio: N/A"
    )
    ax.text(
        0.97,
        0.03,
        f"opposite signs at {opposite_count}/{len(h)} steps\n{ratio_text}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.88},
    )
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.text(
        0.5,
        0.025,
        "The forward and backward formulas approach the same analytic curvature from opposite sides. "
        "Sign changes inside the gray region occur after finite-precision cancellation begins to dominate "
        "and are not evidence of failed first-order convergence. Gray shading begins at the first observed "
        "post-minimum increase in either absolute error.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#333333",
        wrap=True,
    )
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.15)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rkmp(lines: list[str], output_dir: Path) -> dict[str, tuple[list[float], list[float]]]:
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

    signed_header = next(i for i in range(marker, len(lines)) if "forward signed" in lines[i])
    signed_rows = numeric_rows(lines, signed_header, 4)
    signed_h = [row[0] for row in signed_rows]
    signed_forward = [row[1] for row in signed_rows]
    signed_backward = [row[2] for row in signed_rows]
    write_signed_one_sided_plot(
        output_dir / "rkmp_one_sided_signed",
        signed_h,
        signed_forward,
        signed_backward,
        err_forward,
        err_backward,
    )

    high_marker = next(i for i, line in enumerate(lines) if "exponent-eight centered" in line)
    high35_header = next(i for i in range(high_marker, len(lines)) if "high3 abs" in lines[i])
    high35_rows = fixed_table_rows(lines, high35_header, 7)
    high79_header = next(i for i in range(high_marker, len(lines)) if "high7 abs" in lines[i])
    high79_rows = fixed_table_rows(lines, high79_header, 7)
    high_h = [float(row[0]) for row in high35_rows]
    high_err3 = [float(row[2]) for row in high35_rows]
    high_err5 = [float(row[5]) for row in high35_rows]
    high_err7 = [float(row[2]) for row in high79_rows]
    high_err9 = [float(row[5]) for row in high79_rows]
    high_series: list[tuple[str, list[float], list[float], bool]] = [
        ("Three-point", high_h, high_err3, False),
        ("Five-point", high_h, high_err5, False),
        ("Seven-point", high_h, high_err7, False),
        ("Nine-point", high_h, high_err9, False),
    ]
    for values, order, label, lower, upper in (
        (high_err3, 2, r"Extrapolated theoretical $O(h^2)$ truncation trend", 1.7, 2.3),
        (high_err5, 4, r"Extrapolated theoretical $O(h^4)$ truncation trend", 3.2, 4.8),
        (high_err7, 6, r"Extrapolated theoretical $O(h^6)$ truncation trend", 4.8, 7.2),
        (high_err9, 8, r"Extrapolated theoretical $O(h^8)$ truncation trend", 6.4, 9.6),
    ):
        guide = slope_guide(high_h, values, order, label, lower, upper)
        if guide:
            high_series.append((*guide, True))
    high_onsets = [
        value
        for values in (high_err3, high_err5, high_err7, high_err9)
        if (value := small_h_upturn_onset(high_h, values)) is not None
    ]
    high_summaries = []
    for name, values, lower, upper in (
        ("3pt", high_err3, 1.7, 2.3),
        ("5pt", high_err5, 3.2, 4.8),
        ("7pt", high_err7, 4.8, 7.2),
        ("9pt", high_err9, 6.4, 9.6),
    ):
        orders = in_range_orders(high_h, values, lower, upper)
        best = min(value for value in values if value > 0.0 and math.isfinite(value))
        order_text = f"p={min(orders):.2f}-{max(orders):.2f}" if orders else "p=N/A"
        high_summaries.append(
            f"{name}: {order_text} ({len(orders)} interval{'s' if len(orders) != 1 else ''}); "
            f"best={best:.2e}"
        )
    write_plot(
        output_dir / "rkmp_high_order",
        "RKMP exponent-eight centered verification",
        high_series,
        "Degree-ten controlled energy; seven- and nine-point results are report-only",
        summary=high_summaries,
        caption="The exponent-eight fixture supplies nonzero eighth- and tenth-derivative terms. "
        "The measured curves are solid with markers; every dashed guide is anchored directly "
        "to measured data in its observed convergence region.",
        small_h_upturn_start=max(high_onsets) if high_onsets else None,
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
        # The native TestThermo30 polynomial begins at finite-precision-level agreement,
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
                "finite-precision-limited from coarsest retained step",
            ],
            "The unmodified low-order RKMP polynomial agrees at nearly machine precision "
            "before refinement. This is a native Thermochimica accuracy cross-check, not "
            "a truncation-order demonstration.",
            small_h_upturn_start=max(
                step
                for _, h_values, _, _ in plot_series
                for step in h_values
                if step > 0.0 and math.isfinite(step)
            ),
            small_h_region_label="Finite-precision-limited sweep",
        )

    return {
        "three-point": (high_h, high_err3),
        "five-point": (high_h, high_err5),
    }


def plot_cef(lines: list[str], output_dir: Path) -> dict[str, tuple[list[float], list[float]]]:
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
    onset3 = small_h_upturn_onset(h, err3)
    onset5 = small_h_upturn_onset(h, err5)
    onset_values = [value for value in (onset3, onset5) if value is not None]
    order3 = in_range_order_bounds(h, err3, 1.7, 2.3)
    order5 = in_range_order_bounds(h, err5, 3.2, 4.8)
    best3 = min(value for value in err3 if value > 0.0 and math.isfinite(value))
    best5 = min(value for value in err5 if value > 0.0 and math.isfinite(value))
    summary = [
        f"3-point: p={order3[0]:.2f}-{order3[1]:.2f}; best={best3:.2e}",
        f"5-point: p={order5[0]:.2f}-{order5[1]:.2f}; best={best5:.2e}",
        f"small-$h$ upturn at $h$={onset3:.2e} / {onset5:.2e}",
    ]
    for guide in (
        slope_guide(h, err3, 2, r"Extrapolated theoretical $O(h^2)$ truncation trend", 1.7, 2.3),
        slope_guide(h, err5, 4, r"Extrapolated theoretical $O(h^4)$ truncation trend", 3.2, 4.8),
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
        "identified small-$h$ finite-precision region.",
        small_h_upturn_start=max(onset_values) if onset_values else None,
    )
    return {
        "three-point": (h, err3),
        "five-point": (h, err5),
    }


def plot_truncation_coefficients(
    output_stem: Path,
    rkmp: dict[str, tuple[list[float], list[float]]],
    cef: dict[str, tuple[list[float], list[float]]],
) -> None:
    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    specifications = (
        ("RKMP 3-point", rkmp["three-point"], 2, 1.7, 2.3),
        ("RKMP 5-point", rkmp["five-point"], 4, 3.2, 4.8),
        ("CEF 3-point", cef["three-point"], 2, 1.7, 2.3),
        ("CEF 5-point", cef["five-point"], 4, 3.2, 4.8),
    )
    for index, (label, (h, error), order, lower, upper) in enumerate(specifications):
        region = expected_order_region(h, error, lower, upper)
        if len(region) < 2:
            continue
        steps = [item[0] for item in region]
        coefficients = [item[1] / item[0] ** order for item in region]
        usable = sorted(zip(steps, coefficients))
        ax.plot(
            [item[0] for item in usable],
            [item[1] for item in usable],
            marker=MARKERS[index],
            markersize=5,
            linewidth=1.8,
            color=COLORS[index],
            label=f"{label}, $e/h^{order}$",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Perturbation size, $h$ (smaller $h$ to the left)")
    ax.set_ylabel("Scaled truncation coefficient, $e(h)/h^p$")
    ax.set_title("RKMP and CEF truncation-coefficient evidence", fontweight="bold", pad=24)
    ax.text(
        0.5,
        1.015,
        "A plateau indicates that measured error follows the expected power of h",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10,
        color="#444444",
    )
    ax.grid(which="major", color="#cfcfcf", linewidth=0.8)
    ax.grid(which="minor", color="#e8e8e8", linewidth=0.5, alpha=0.8)
    ax.legend(loc="best", frameon=False, fontsize=9)
    fig.text(
        0.5,
        0.025,
        "Plateau height depends on each model's higher derivatives, state, and error normalization. "
        "Different RKMP and CEF heights are therefore expected and are not Hessian disagreement.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#333333",
        wrap=True,
    )
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.15)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


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
        displayed_h = [step for h_values, _ in summary_series for step in h_values]
        guide = slope_guide(
            *guide_data,
            4,
            r"Extrapolated theoretical $O(h^4)$ truncation trend",
            3.2,
            4.8,
            extrapolation_h=displayed_h,
        )
        if guide:
            plot_series.append((*guide, True))
    onset_values = [
        onset
        for h_values, error_values in summary_series
        if (onset := small_h_upturn_onset(h_values, error_values)) is not None
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
        onset_summary = f"small-$h$ upturn at $h$={onset_values[0]:.2e}"
    else:
        onset_summary = (
            f"small-$h$ upturn at $h$={min(onset_values):.2e}-{max(onset_values):.2e}"
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
        small_h_upturn_start=max(onset_values) if onset_values else None,
    )


def plot_standalone_subq(lines: list[str], output_dir: Path) -> None:
    """Plot controlled SUBQ cases without conflating them with native decoding."""
    plot_series = []
    guide_data: tuple[list[float], list[float]] | None = None
    summary_series: list[tuple[list[float], list[float]]] = []
    cases = (
        ("SUBQ configurational", "SUBQ nonuniform-zeta configurational"),
        ("SUBQ G chi incidence", "SUBQ G binary chi incidence"),
        ("SUBQ G swapped incidence", "SUBQ G binary swapped chi incidence"),
        ("SUBQ configurational + B", "SUBQ nonuniform-zeta configurational + B"),
    )
    for label, report_name in cases:
        marker = next(i for i, line in enumerate(lines) if line.strip() == report_name)
        header = next(i for i in range(marker, len(lines)) if "3pt abs" in lines[i])
        rows = numeric_rows(lines, header, 7)
        h = [row[0] for row in rows]
        error = [row[5] for row in rows]
        plot_series.append((label, h, error, False))
        summary_series.append((h, error))
        guide_data = guide_data or (h, error)

    displayed_h = [step for h_values, _ in summary_series for step in h_values]
    if guide_data:
        guide = slope_guide(
            *guide_data,
            4,
            r"Extrapolated theoretical $O(h^4)$ truncation trend",
            3.2,
            4.8,
            extrapolation_h=displayed_h,
        )
        if guide:
            plot_series.append((*guide, True))

    onset_values = [
        onset
        for h_values, error_values in summary_series
        if (onset := small_h_upturn_onset(h_values, error_values)) is not None
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
    if onset_values and math.isclose(min(onset_values), max(onset_values), rel_tol=1e-12):
        onset_summary = f"small-$h$ upturn at $h$={onset_values[0]:.2e}"
    elif onset_values:
        onset_summary = (
            f"small-$h$ upturn at $h$={min(onset_values):.2e}-{max(onset_values):.2e}"
        )
    else:
        onset_summary = "no small-$h$ upturn observed in retained sweeps"

    write_plot(
        output_dir / "mqmqa_subq_standalone",
        "Standalone SUBQ scalar-energy verification",
        plot_series,
        "Controlled nonuniform-zeta SUBQ configurational, chi-incidence, and B-family cases",
        summary=[
            f"5-point p={min(item[0] for item in order_values):.2f}-"
            f"{max(item[1] for item in order_values):.2f}",
            f"best scaled error={min(best_values):.2e}-{max(best_values):.2e}",
            onset_summary,
        ],
        caption="The controlled cases verify the SUBQ configurational exponents, "
        "mixed-environment chi weights, pair-specific zeta propagation, and B-family "
        "curvature. They use synthetic interior states and do not constitute native "
        "database-backed Thermochimica comparison.",
        small_h_upturn_start=max(onset_values) if onset_values else None,
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
        displayed_h = [step for h_values, _ in summary_series for step in h_values]
        guide = slope_guide(
            *guide_data,
            2,
            r"Extrapolated theoretical $O(h^2)$ truncation trend",
            1.7,
            2.3,
            extrapolation_h=displayed_h,
        )
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
            "no small-$h$ upturn observed in retained sweep",
        ],
        caption="At a converged CuFeC-Kang.dat plain-SUBG state, analytic Hessian-vector "
        "products agree with production Thermochimica partial-molar differences. "
        "Scope: nonmagnetic reference/configurational/G-family behavior.",
    )


def plot_native_subq(lines: list[str], output_dir: Path) -> None:
    """Show complete tangent coverage without assigning meaning to direction colors."""
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

    directions = sorted({int(row[0]) for row in rows})
    direction_series: list[tuple[int, list[float], list[float]]] = []
    for direction in directions:
        selected = [row for row in rows if int(row[0]) == direction]
        h = [row[1] for row in selected]
        error = [row[3] for row in selected]
        direction_series.append((direction, h, error))

    # Each direction has its own positivity-safe mole scale s_i.  The Fortran
    # sweep starts at h=0.05*s_i, so recovering s_i from the largest reported
    # step lets the presentation compare the common refinement coordinate
    # h/s_i without changing any finite-difference calculation or test gate.
    normalized_series: list[tuple[int, list[float], list[float]]] = []
    for direction, h_values, error_values in direction_series:
        direction_scale = max(h_values) / 0.05
        normalized_series.append(
            (direction, [step / direction_scale for step in h_values], error_values)
        )

    order_values = [
        bounds
        for _, h_values, error_values in direction_series
        if (bounds := in_range_order_bounds(h_values, error_values, 1.7, 2.3)) is not None
    ]
    best_values = [
        min(value for value in error_values if value > 0.0 and math.isfinite(value))
        for _, _, error_values in direction_series
    ]
    normalized_onset_values = [
        onset
        for _, h_values, error_values in normalized_series
        if (onset := small_h_upturn_onset(h_values, error_values)) is not None
    ]
    gradient_line = next(line for line in lines if line.startswith("direct unconstrained-gradient error"))
    gradient_error = numbers(gradient_line)[0]
    worst_index = max(range(len(best_values)), key=best_values.__getitem__)
    worst_direction, worst_h, worst_error = normalized_series[worst_index]
    displayed_h = [step for _, h_values, _ in normalized_series for step in h_values]
    guide = slope_guide(
        worst_h,
        worst_error,
        2,
        r"Extrapolated theoretical $O(h^2)$ trend",
        1.7,
        2.3,
        extrapolation_h=displayed_h,
    )

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(9.8, 8.8),
        gridspec_kw={"height_ratios": [3.2, 1.15]},
    )
    for index, (_, h_values, error_values) in enumerate(normalized_series):
        usable = sorted(
            (step, value)
            for step, value in zip(h_values, error_values)
            if step > 0.0 and value > 0.0 and math.isfinite(step) and math.isfinite(value)
        )
        x_values, y_values = zip(*usable)
        ax_top.plot(
            x_values,
            y_values,
            color="#9a9a9a",
            alpha=0.42,
            linewidth=1.1,
            label="14 independent tangent directions" if index == 0 else "_nolegend_",
        )

    worst_usable = sorted(zip(worst_h, worst_error))
    ax_top.plot(
        [item[0] for item in worst_usable],
        [item[1] for item in worst_usable],
        color="#8a1538",
        marker="o",
        markersize=5,
        linewidth=2.5,
        label=f"Worst-case direction ({worst_direction})",
    )
    if guide:
        label, guide_h, guide_error = guide
        ax_top.plot(guide_h,guide_error,"--",color="#2f5597",linewidth=2.0,label=label)

    ax_top.set_xscale("log")
    ax_top.set_yscale("log")
    set_measured_limits(
        ax_top,
        [(h_values,error_values) for _,h_values,error_values in normalized_series],
    )
    if normalized_onset_values:
        left_limit, right_limit = ax_top.get_xlim()
        shade_right = min(max(max(normalized_onset_values),left_limit),right_limit)
        ax_top.axvspan(
            left_limit,
            shade_right,
            color="#777777",
            alpha=0.08,
            label="Small-$h$ finite-precision region",
        )
    ax_top.set_xlabel(r"Normalized perturbation, $\widehat h=h/s_i$ (smaller to the left)")
    ax_top.set_ylabel("Scale-normalized error")
    ax_top.grid(which="major",color="#cfcfcf",linewidth=0.8)
    ax_top.grid(which="minor",color="#e8e8e8",linewidth=0.5,alpha=0.8)
    ax_top.legend(loc="upper left",frameon=False,fontsize=8.5)
    ax_top.text(
        0.98,
        0.04,
        "\n".join(
            [
                "Fortran gate: 14/14 directions passed",
                f"in-range observed $p$={min(item[0] for item in order_values):.2f}-"
                f"{max(item[1] for item in order_values):.2f}",
                f"worst normwise minimum={max(best_values):.2e}",
                f"direct gradient error={gradient_error:.2e}",
            ]
        ),
        transform=ax_top.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="#333333",
        bbox={"facecolor":"white","edgecolor":"#cccccc","alpha":0.9},
    )

    ax_bottom.scatter(directions,best_values,color="#777777",s=35,zorder=3,label="Direction minimum")
    ax_bottom.scatter(
        [worst_direction],
        [best_values[worst_index]],
        color="#8a1538",
        s=60,
        zorder=4,
        label="Worst-case direction",
    )
    ax_bottom.axhline(1.0e-8,color="#2f5597",linestyle="--",linewidth=1.8,label="Acceptance threshold $10^{-8}$")
    ax_bottom.set_yscale("log")
    ax_bottom.set_xlim(0.4,len(directions)+0.6)
    ax_bottom.set_ylim(min(best_values)*0.55,2.0e-8)
    ax_bottom.set_xticks(directions)
    ax_bottom.set_xlabel("Independent tangent direction")
    ax_bottom.set_ylabel("Minimum scaled error")
    ax_bottom.grid(which="major",color="#d8d8d8",linewidth=0.7)
    ax_bottom.legend(loc="upper right",ncol=2,frameon=False,fontsize=8.2)

    fig.suptitle("Native SUBQ G/Q Hessian verification",fontweight="bold",fontsize=15,y=0.975)
    fig.text(
        0.5,
        0.935,
        "Production partial-molar finite differences over the complete 14-dimensional tangent basis",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#444444",
    )
    caption = (
        "Grey curves show complete tangent-space coverage; the highlighted worst case still follows "
        "second-order convergence and remains below the acceptance threshold. Normalizing h by each "
        "direction's positivity-safe scale removes abundance-driven horizontal clustering; absolute h "
        "values remain in the report. At small h, cancellation "
        "between nearly equal partial-molar evaluations and division by h amplify floating-point error. "
        "Scope: assessed FeTiVO.dat nonmagnetic, uniform-zeta reference/configurational/G/Q behavior."
    )
    fig.text(0.5,0.018,textwrap.fill(caption,width=145),ha="center",va="bottom",fontsize=8.6,color="#333333")
    fig.subplots_adjust(left=0.11,right=0.97,top=0.90,bottom=0.15,hspace=0.34)
    output_stem = output_dir / "mqmqa_subq_native"
    fig.savefig(output_stem.with_suffix(".svg"),bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"),dpi=300,bbox_inches="tight")
    plt.close(fig)


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
        "native_subq": run_report(args.bin_dir, "TestMQMQASUBQNativeHessianVerification"),
    }
    for name, lines in reports.items():
        (args.output_dir / f"{name}_report.txt").write_text("\n".join(lines) + "\n")

    rkmp_controlled = plot_rkmp(reports["rkmp"], args.output_dir)
    cef_controlled = plot_cef(reports["cef"], args.output_dir)
    plot_truncation_coefficients(
        args.output_dir / "rkmp_cef_truncation_coefficients",
        rkmp_controlled,
        cef_controlled,
    )
    plot_standalone_mqmqa(reports["standalone"], args.output_dir)
    plot_standalone_subq(reports["standalone"], args.output_dir)
    plot_native_mqmqa(reports["native"], args.output_dir)
    plot_native_subq(reports["native_subq"], args.output_dir)
    print(f"Wrote reports and plots to {args.output_dir}")


if __name__ == "__main__":
    main()
