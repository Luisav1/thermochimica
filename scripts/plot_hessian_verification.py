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
from html import escape
from pathlib import Path

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
COLORS = ("#8a1538", "#006747", "#2f5597", "#c55a11", "#7030a0", "#444444")


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
    guide_e = [e0 * (value / h0) ** order for value in guide_h]
    return label, guide_h, guide_e


def write_svg_plot(
    output: Path,
    title: str,
    series: list[tuple[str, list[float], list[float], bool]],
) -> None:
    width, height = 900, 560
    left, right, top, bottom = 90, 30, 55, 105
    plot_width = width - left - right
    plot_height = height - top - bottom
    points = [
        (x, y)
        for _, xs, ys, _ in series
        for x, y in zip(xs, ys)
        if x > 0.0 and y > 0.0 and math.isfinite(x) and math.isfinite(y)
    ]
    xmin, xmax = min(math.log10(x) for x, _ in points), max(math.log10(x) for x, _ in points)
    ymin, ymax = min(math.log10(y) for _, y in points), max(math.log10(y) for _, y in points)
    xpad = max(0.15, 0.04 * (xmax - xmin))
    ypad = max(0.25, 0.04 * (ymax - ymin))
    xmin, xmax, ymin, ymax = xmin - xpad, xmax + xpad, ymin - ypad, ymax + ypad

    def sx(value: float) -> float:
        return left + (math.log10(value) - xmin) * plot_width / (xmax - xmin)

    def sy(value: float) -> float:
        return top + (ymax - math.log10(value)) * plot_height / (ymax - ymin)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="Arial" '
        f'font-size="20" font-weight="bold">{escape(title)}</text>',
    ]
    for exponent in range(math.floor(xmin), math.ceil(xmax) + 1):
        x = sx(10.0**exponent)
        svg.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top+plot_height}" '
                   'stroke="#dddddd" stroke-width="1"/>')
        svg.append(f'<text x="{x:.1f}" y="{top+plot_height+22}" text-anchor="middle" '
                   f'font-family="Arial" font-size="12">10^{exponent}</text>')
    for exponent in range(math.floor(ymin), math.ceil(ymax) + 1):
        y = sy(10.0**exponent)
        svg.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_width}" y2="{y:.1f}" '
                   'stroke="#dddddd" stroke-width="1"/>')
        svg.append(f'<text x="{left-10}" y="{y+4:.1f}" text-anchor="end" '
                   f'font-family="Arial" font-size="12">10^{exponent}</text>')
    svg.extend([
        f'<line x1="{left}" y1="{top+plot_height}" x2="{left+plot_width}" '
        f'y2="{top+plot_height}" stroke="black" stroke-width="1.5"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_height}" '
        'stroke="black" stroke-width="1.5"/>',
        f'<text x="{left+plot_width/2}" y="{height-56}" text-anchor="middle" '
        'font-family="Arial" font-size="14">Perturbation size, h</text>',
        f'<text x="22" y="{top+plot_height/2}" text-anchor="middle" font-family="Arial" '
        'font-size="14" transform="rotate(-90 22 '
        f'{top+plot_height/2})">Scale-normalized error</text>',
    ])
    legend_y = height - 31
    legend_x = left
    for index, (label, xs, ys, dashed) in enumerate(series):
        color = COLORS[index % len(COLORS)]
        usable = [(sx(x), sy(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
        if len(usable) < 2:
            continue
        path = " ".join(f"{'M' if i == 0 else 'L'} {x:.2f} {y:.2f}" for i, (x, y) in enumerate(usable))
        dash = ' stroke-dasharray="7,5"' if dashed else ""
        svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2"{dash}/>')
        if not dashed:
            for x, y in usable:
                svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}"/>')
        item_x = legend_x + (index % 3) * 250
        item_y = legend_y + (index // 3) * 22
        svg.append(f'<line x1="{item_x}" y1="{item_y}" x2="{item_x+25}" y2="{item_y}" '
                   f'stroke="{color}" stroke-width="2"{dash}/>')
        svg.append(f'<text x="{item_x+31}" y="{item_y+4}" font-family="Arial" '
                   f'font-size="11">{escape(label)}</text>')
    svg.append("</svg>")
    output.write_text("\n".join(svg) + "\n")


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
    for guide in (slope_guide(h, err3, 2, "O(h^2) guide"), slope_guide(h, err5, 4, "O(h^4) guide")):
        if guide:
            plot_series.append((*guide, True))
    write_svg_plot(output_dir / "rkmp_energy.svg", "RKMP controlled excess-energy verification", plot_series)

    plot_series = []
    cursor = 0
    direction = 0
    while True:
        try:
            marker, rows = find_after(lines, "direction: species", "norm abs", cursor)
        except StopIteration:
            break
        direction += 1
        cursor = marker + 1
        h = [row[0] for row in rows]
        error = [row[2] for row in rows]
        if any(value > 0.0 for value in error):
            plot_series.append((f"production partial-molar direction {direction}", h, error, False))
    if plot_series:
        guide = slope_guide(plot_series[0][1], plot_series[0][2], 2, "O(h^2) guide")
        if guide:
            plot_series.append((*guide, True))
        write_svg_plot(output_dir / "rkmp_partial_molar.svg", "RKMP production partial-molar verification", plot_series)


def plot_cef(lines: list[str], output_dir: Path) -> None:
    marker = next(i for i, line in enumerate(lines) if line.strip() == "ALABANDITE controlled")
    header = next(i for i in range(marker, len(lines)) if "3pt abs" in lines[i])
    rows = fixed_table_rows(lines, header, 8)
    rows = [row for row in rows if int(row[0]) == 1]
    h = [float(row[1]) for row in rows]
    err3 = [float(row[3]) for row in rows]
    err5 = [float(row[6]) for row in rows]
    plot_series = [
        ("CEF ALABANDITE controlled, 3-point", h, err3, False),
        ("CEF ALABANDITE controlled, 5-point", h, err5, False),
    ]
    for guide in (slope_guide(h, err3, 2, "O(h^2) guide"), slope_guide(h, err5, 4, "O(h^4) guide")):
        if guide:
            plot_series.append((*guide, True))
    write_svg_plot(output_dir / "cef_energy.svg", "CEF controlled scalar-energy verification", plot_series)


def plot_standalone_mqmqa(lines: list[str], output_dir: Path) -> None:
    plot_series = []
    guide_data: tuple[list[float], list[float]] | None = None
    for case in ("G binary", "Q binary", "B family"):
        marker = next(i for i, line in enumerate(lines) if line.strip() == case)
        header = next(i for i in range(marker, len(lines)) if "3pt abs" in lines[i])
        rows = numeric_rows(lines, header, 7)
        h = [row[0] for row in rows]
        error = [row[5] for row in rows]
        plot_series.append((f"standalone MQMQA {case}, 5-point", h, error, False))
        guide_data = guide_data or (h, error)
    if guide_data:
        guide = slope_guide(*guide_data, 4, "O(h^4) guide")
        if guide:
            plot_series.append((*guide, True))
    write_svg_plot(output_dir / "mqmqa_standalone.svg", "Standalone MQMQA G/Q/B verification", plot_series)


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
    for direction in sorted({int(row[0]) for row in rows}):
        selected = [row for row in rows if int(row[0]) == direction]
        h = [row[1] for row in selected]
        error = [row[3] for row in selected]
        plot_series.append((f"native MQMQA direction {direction}", h, error, False))
        guide_data = guide_data or (h, error)
    if guide_data:
        guide = slope_guide(*guide_data, 2, "O(h^2) guide")
        if guide:
            plot_series.append((*guide, True))
    write_svg_plot(output_dir / "mqmqa_native.svg", "Native MQMQA G-family verification", plot_series)


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
