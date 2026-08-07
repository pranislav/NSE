#!/usr/bin/env python3
"""
Visualize MMS solutions and finer integral error fields from legacy ASCII VTK files.

The script expects folders named like

    solns/mms_re100_deg1_finer_integral_error
    solns/mms_re100_deg2_finer_integral_error
    solns/mms_re100_deg3_finer_integral_error

and, when present, the analogous Re 7500 folders.  For each folder it selects
the VTK file with the highest number after "ref"; ties are resolved by the
highest number after "newt".
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np

from visualize_mms_errors import (
    DEGREES,
    VtkData,
    field_values,
    make_triangles,
    plot_panel,
    read_legacy_ascii_vtk,
    select_vtk_file,
)


@dataclass(frozen=True)
class IntegralErrorRow:
    error_array: str
    label: str
    solution_field: str
    re_value: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build an MMS solution/finer-integral-error visualization from VTK files."
        )
    )
    parser.add_argument("--solns-dir", type=Path, default=Path("solns"))
    parser.add_argument(
        "--output", type=Path, default=Path("mms_integral_errors.png")
    )
    parser.add_argument(
        "--folder-suffix",
        default="_finer_integral_error",
        help="Suffix appended to the old mms_re*_deg* directory names.",
    )
    parser.add_argument(
        "--re-values",
        type=int,
        nargs="+",
        default=None,
        help=(
            "Re values to plot. Defaults to all Re values with complete deg1-deg3 "
            "folders for the selected suffix."
        ),
    )
    parser.add_argument(
        "--solution-degree",
        type=int,
        default=3,
        choices=DEGREES,
        help="Degree folder used for the original solution column.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output DPI. Default gives a large image suitable for zooming.",
    )
    parser.add_argument("--show", action="store_true", help="Show the figure interactively.")
    return parser.parse_args()


def folder_for(solns_dir: Path, re_value: int, degree: int, suffix: str) -> Path:
    return solns_dir / f"mms_re{re_value}_deg{degree}{suffix}"


def discover_re_values(solns_dir: Path, suffix: str) -> list[int]:
    available: dict[int, set[int]] = {}
    for folder in solns_dir.iterdir():
        if not folder.is_dir() or not folder.name.endswith(suffix):
            continue
        prefix = folder.name[: -len(suffix)] if suffix else folder.name
        match = re.fullmatch(r"mms_re(?P<re>\d+)_deg(?P<degree>\d+)", prefix)
        if match is None:
            continue
        re_value = int(match.group("re"))
        degree = int(match.group("degree"))
        available.setdefault(re_value, set()).add(degree)

    return sorted(
        re_value
        for re_value, degrees in available.items()
        if all(degree in degrees for degree in DEGREES)
    )


def exact_array_values(data: VtkData, name: str) -> np.ndarray:
    if name not in data.point_arrays:
        raise KeyError(f"Array '{name}' not found; available: {sorted(data.point_arrays)}")

    values = data.point_arrays[name]
    if values.ndim == 2:
        return np.linalg.norm(values[:, :2], axis=1)
    return values


def row_specs(re_values: list[int]) -> list[IntegralErrorRow]:
    field_specs = [
        ("velocity_L2_cell_error", "velocity L2", "velocity"),
        ("velocity_H1_cell_error", "velocity H1", "velocity"),
        ("pressure_L2_cell_error", "pressure L2", "pressure"),
    ]
    return [
        IntegralErrorRow(
            error_array=error_array,
            label=label,
            solution_field=solution_field,
            re_value=re_value,
        )
        for error_array, label, solution_field in field_specs
        for re_value in re_values
    ]


def main() -> None:
    args = parse_args()

    re_values = args.re_values or discover_re_values(args.solns_dir, args.folder_suffix)
    if not re_values:
        raise FileNotFoundError(
            f"No complete deg1-deg3 MMS folders found in {args.solns_dir} "
            f"with suffix '{args.folder_suffix}'"
        )

    selected: dict[tuple[int, int], Path] = {}
    vtk_data: dict[tuple[int, int], VtkData] = {}
    triangles: dict[tuple[int, int], np.ndarray] = {}

    for re_value in re_values:
        for degree in DEGREES:
            folder = folder_for(args.solns_dir, re_value, degree, args.folder_suffix)
            path = select_vtk_file(folder)
            data = read_legacy_ascii_vtk(path)
            selected[(re_value, degree)] = path
            vtk_data[(re_value, degree)] = data
            triangles[(re_value, degree)] = make_triangles(data.cells)

    columns = [
        f"solution deg{args.solution_degree}",
        "error deg1",
        "error deg2",
        "error deg3",
    ]
    rows = row_specs(re_values)

    row_count = len(rows)
    fig, axes = plt.subplots(row_count, 4, figsize=(28, 5.5 * row_count), constrained_layout=False)
    axes = np.atleast_2d(axes)
    fig.subplots_adjust(left=0.1, right=0.985, bottom=0.035, top=0.88, wspace=0.18, hspace=0.08)
    fig.suptitle("MMS solutions and cell (integral) errors", fontsize=34, y=0.975)

    for row, spec in enumerate(rows):
        for col in range(4):
            degree = args.solution_degree if col == 0 else col
            data = vtk_data[(spec.re_value, degree)]
            values = (
                field_values(data, spec.solution_field, error=False)
                if col == 0
                else exact_array_values(data, spec.error_array)
            )
            plot_panel(fig, axes[row, col], data, triangles[(spec.re_value, degree)], values)

    fig.canvas.draw()
    for col, title in enumerate(columns):
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

    for row, spec in enumerate(rows):
        bbox = axes[row, 0].get_position()
        fig.text(
            0.045,
            (bbox.y0 + bbox.y1) / 2,
            f"{spec.label}\nRe {spec.re_value}",
            ha="center",
            va="center",
            fontsize=23,
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
