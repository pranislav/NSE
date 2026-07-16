#!/usr/bin/env python3
"""
Visualize MMS pressure/velocity solutions and errors from legacy ASCII VTK files.

The script expects folders named like

    solns/mms_re100_deg1
    solns/mms_re100_deg2
    solns/mms_re100_deg3
    solns/mms_re7500_deg1
    solns/mms_re7500_deg2
    solns/mms_re7500_deg3

For each folder it selects the VTK file with the highest number after "ref";
ties are resolved by the highest number after "newt".
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import re
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


RE_VALUES = (100, 7500)
DEGREES = (1, 2, 3)
REF_NEWT_RE = re.compile(r"ref(?P<ref>\d+).*newt(?P<newt>\d+)", re.IGNORECASE)


@dataclass(frozen=True)
class VtkData:
    points: np.ndarray
    cells: list[list[int]]
    point_arrays: dict[str, np.ndarray]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a 4x4 MMS solution/error visualization from VTK files."
    )
    parser.add_argument("--solns-dir", type=Path, default=Path("solns"))
    parser.add_argument("--output", type=Path, default=Path("mms_errors.png"))
    parser.add_argument(
        "--solution-degree",
        type=int,
        default=3,
        choices=DEGREES,
        help="Degree folder used for the solution column.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output DPI. Default gives a large image suitable for zooming.",
    )
    parser.add_argument("--show", action="store_true", help="Show the figure interactively.")
    return parser.parse_args()


def vtk_sort_key(path: Path) -> tuple[int, int, str] | None:
    match = REF_NEWT_RE.search(path.name)
    if match is None:
        return None
    return (int(match.group("ref")), int(match.group("newt")), path.name)


def select_vtk_file(folder: Path) -> Path:
    candidates: list[tuple[tuple[int, int, str], Path]] = []
    for path in folder.glob("*.vtk"):
        key = vtk_sort_key(path)
        if key is not None:
            candidates.append((key, path))

    if not candidates:
        raise FileNotFoundError(f"No ref*/newt* VTK files found in {folder}")

    return max(candidates, key=lambda item: item[0])[1]


def read_numbers(lines: list[str], start: int, count: int) -> tuple[np.ndarray, int]:
    chunks: list[str] = []
    got = 0
    index = start
    while got < count and index < len(lines):
        line = lines[index].strip()
        if line:
            chunks.append(line)
            got += len(line.split())
        index += 1

    values = np.fromstring(" ".join(chunks), sep=" ")
    if values.size < count:
        raise ValueError(f"Expected {count} numeric values, found {values.size}")

    return values[:count], index


def read_legacy_ascii_vtk(path: Path) -> VtkData:
    lines = path.read_text().splitlines()
    points: np.ndarray | None = None
    cells: list[list[int]] = []
    point_arrays: dict[str, np.ndarray] = {}
    point_count: int | None = None

    index = 0
    while index < len(lines):
        words = lines[index].strip().split()
        if not words:
            index += 1
            continue

        keyword = words[0].upper()
        if keyword == "POINTS":
            point_count = int(words[1])
            values, index = read_numbers(lines, index + 1, point_count * 3)
            points = values.reshape(point_count, 3)[:, :2]
            continue

        if keyword == "CELLS":
            cell_count = int(words[1])
            cells = []
            index += 1
            for _ in range(cell_count):
                cell_words = lines[index].strip().split()
                vertex_count = int(cell_words[0])
                cells.append([int(value) for value in cell_words[1 : 1 + vertex_count]])
                index += 1
            continue

        if keyword == "POINT_DATA":
            point_count = int(words[1])
            index += 1
            continue

        if keyword == "VECTORS":
            if point_count is None:
                raise ValueError(f"VECTORS before POINT_DATA in {path}")
            name = words[1]
            values, index = read_numbers(lines, index + 1, point_count * 3)
            point_arrays[name] = values.reshape(point_count, 3)
            continue

        if keyword == "SCALARS":
            if point_count is None:
                raise ValueError(f"SCALARS before POINT_DATA in {path}")
            name = words[1]
            components = int(words[3]) if len(words) > 3 else 1
            index += 1
            if index < len(lines) and lines[index].strip().upper().startswith("LOOKUP_TABLE"):
                index += 1
            values, index = read_numbers(lines, index, point_count * components)
            array = values.reshape(point_count, components)
            point_arrays[name] = array[:, 0] if components == 1 else array
            continue

        index += 1

    if points is None:
        raise ValueError(f"No POINTS section found in {path}")
    if not cells:
        raise ValueError(f"No CELLS section found in {path}")

    return VtkData(points=points, cells=cells, point_arrays=point_arrays)


def make_triangles(cells: Iterable[Iterable[int]]) -> np.ndarray:
    triangles: list[list[int]] = []
    for cell in cells:
        ids = list(cell)
        if len(ids) == 3:
            triangles.append(ids)
        elif len(ids) == 4:
            triangles.append([ids[0], ids[1], ids[2]])
            triangles.append([ids[0], ids[2], ids[3]])
        elif len(ids) > 4:
            for offset in range(1, len(ids) - 1):
                triangles.append([ids[0], ids[offset], ids[offset + 1]])
    return np.asarray(triangles, dtype=int)


def field_values(data: VtkData, field: str, error: bool) -> np.ndarray:
    name = f"{field}_error" if error else field
    if name not in data.point_arrays:
        raise KeyError(f"Array '{name}' not found; available: {sorted(data.point_arrays)}")

    values = data.point_arrays[name]
    if values.ndim == 2:
        return np.linalg.norm(values[:, :2], axis=1)
    return values


def plot_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    data: VtkData,
    triangles: np.ndarray,
    values: np.ndarray,
) -> None:
    triangulation = mtri.Triangulation(data.points[:, 0], data.points[:, 1], triangles)
    image = ax.tripcolor(triangulation, values, shading="gouraud", cmap="coolwarm")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(np.min(data.points[:, 0])), float(np.max(data.points[:, 0])))
    ax.set_ylim(float(np.min(data.points[:, 1])), float(np.max(data.points[:, 1])))
    ax.margins(0.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.01)
    colorbar.ax.tick_params(labelsize=11)


def folder_for(solns_dir: Path, re_value: int, degree: int) -> Path:
    return solns_dir / f"mms_re{re_value}_deg{degree}"


def main() -> None:
    args = parse_args()

    selected: dict[tuple[int, int], Path] = {}
    vtk_data: dict[tuple[int, int], VtkData] = {}
    triangles: dict[tuple[int, int], np.ndarray] = {}

    for re_value in RE_VALUES:
        for degree in DEGREES:
            folder = folder_for(args.solns_dir, re_value, degree)
            path = select_vtk_file(folder)
            data = read_legacy_ascii_vtk(path)
            selected[(re_value, degree)] = path
            vtk_data[(re_value, degree)] = data
            triangles[(re_value, degree)] = make_triangles(data.cells)

    column_titles = [f"solution deg{args.solution_degree}", "error deg1", "error deg2", "error deg3"]
    row_specs = [
        ("velocity", 100),
        ("velocity", 7500),
        ("pressure", 100),
        ("pressure", 7500),
    ]

    fig, axes = plt.subplots(4, 4, figsize=(28, 22), constrained_layout=False)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.035, top=0.88, wspace=0.18, hspace=0.08)
    fig.suptitle("MMS solutions and errors", fontsize=34, y=0.975)

    for row, (field, re_value) in enumerate(row_specs):
        for col in range(4):
            degree = args.solution_degree if col == 0 else col
            data = vtk_data[(re_value, degree)]
            values = field_values(data, field, error=col != 0)
            plot_panel(fig, axes[row, col], data, triangles[(re_value, degree)], values)

    fig.canvas.draw()
    for col, title in enumerate(column_titles):
        bbox = axes[0, col].get_position()
        fig.text(
            (bbox.x0 + bbox.x1) / 2,
            0.905,
            title,
            ha="center",
            va="bottom",
            fontsize=24,
            fontweight="bold",
        )

    for row, (field, re_value) in enumerate(row_specs):
        bbox = axes[row, 0].get_position()
        fig.text(
            0.035,
            (bbox.y0 + bbox.y1) / 2,
            f"{field}\nRe {re_value}",
            ha="center",
            va="center",
            rotation=90,
            fontsize=24,
            fontweight="bold",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)
    print(f"saved {args.output}")
    print("selected VTK files:")
    for (re_value, degree), path in sorted(selected.items()):
        print(f"  Re {re_value}, deg{degree}: {path}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
