#!/usr/bin/env python3
"""Create a 10-point ternary overlay for the FLiBe-UF4 phase diagram."""

from __future__ import annotations

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = ROOT / "inputs" / "6NE3-individual_project_verify_CaseA.ti"
BASE_IMAGE = ROOT / "flibe-Uf4 ternary phase diagram.png"
OUTPUT_IMAGE = ROOT / "flibe-Uf4 ternary phase diagram 10-point validation overlay.png"

# Pixel coordinates of the ternary vertices in the raster diagram.
# Li at lower left, Be at lower right, U at the apex.
LI_VERTEX = (27.0, 606.0)
BE_VERTEX = (774.0, 622.0)
U_VERTEX = (397.0, 3.0)


def parse_validated_points(path: Path) -> list[tuple[float, float, float, int]]:
    """Return (Li, Be, U, T_liq_C) for the 10 bracket-midpoint compositions."""

    target_temps_c = [440, 430, 480, 580, 625, 675, 770, 827.5, 860, 950]

    lines = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            raw = raw.strip()
            if not raw or raw.startswith("!"):
                continue
            if raw.startswith("nCalc") or raw.startswith("nEl") or raw.startswith("iEl"):
                continue
            parts = raw.split()
            # Data lines are: T[K] printmode Li Be F U
            if len(parts) == 6:
                temp_k = float(parts[0])
                li = float(parts[2])
                be = float(parts[3])
                u = float(parts[5])
                lines.append((temp_k, li, be, u))

    if len(lines) != 20:
        raise ValueError(f"Expected 20 calculation lines, found {len(lines)}")
    if len(target_temps_c) != 10:
        raise ValueError("Expected 10 target temperatures")

    points = []
    for i in range(0, len(lines), 2):
        t1, li1, be1, u1 = lines[i]
        t2, li2, be2, u2 = lines[i + 1]
        if (li1, be1, u1) != (li2, be2, u2):
            raise ValueError(f"Bracket pair {i//2 + 1} does not have matching composition")
        points.append((li1, be1, u1, target_temps_c[i // 2]))

    return points


def barycentric_to_pixel(li: float, be: float, u: float) -> tuple[float, float]:
    total = li + be + u
    f_li = li / total
    f_be = be / total
    f_u = u / total
    x = f_li * LI_VERTEX[0] + f_be * BE_VERTEX[0] + f_u * U_VERTEX[0]
    y = f_li * LI_VERTEX[1] + f_be * BE_VERTEX[1] + f_u * U_VERTEX[1]
    return x, y


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def main() -> None:
    points = parse_validated_points(INPUT_FILE)

    image = Image.open(BASE_IMAGE).convert("RGBA")
    draw = ImageDraw.Draw(image, "RGBA")
    font = load_font(15)

    marker_fill = (255, 102, 0, 255)
    marker_outline = (255, 255, 255, 255)
    label_fill = (28, 62, 124, 255)

    # Small offsets keep labels readable without overwhelming the diagram.
    offsets = [
        (12, -16),
        (12, -16),
        (12, -16),
        (12, -16),
        (12, -16),
        (12, -16),
        (12, -16),
        (-42, -16),
        (-42, -16),
        (-42, -16),
    ]

    for (li, be, u, temp_c), (dx, dy) in zip(points, offsets):
        x, y = barycentric_to_pixel(li, be, u)
        r = 6
        draw.ellipse((x - r, y - r, x + r, y + r), fill=marker_fill, outline=marker_outline, width=2)

        label = f"{temp_c:g} C"
        left, top, right, bottom = draw.textbbox((0, 0), label, font=font, stroke_width=2)
        text_w = right - left
        text_h = bottom - top
        tx = x + dx
        ty = y + dy - text_h / 2

        # White backing keeps the label legible over contour lines.
        pad_x = 3
        pad_y = 2
        draw.rounded_rectangle(
            (tx - pad_x, ty - pad_y, tx + text_w + pad_x, ty + text_h + pad_y),
            radius=3,
            fill=(255, 255, 255, 220),
            outline=(255, 255, 255, 255),
        )
        draw.text(
            (tx, ty),
            label,
            font=font,
            fill=label_fill,
            stroke_width=2,
            stroke_fill=(255, 255, 255, 255),
        )

    image.save(OUTPUT_IMAGE)
    print(OUTPUT_IMAGE)


if __name__ == "__main__":
    main()
