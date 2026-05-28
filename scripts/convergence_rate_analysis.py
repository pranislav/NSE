#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_org_table(path: Path) -> pd.DataFrame:
    """
    Parse dealii convergence table exported in org format.
    """

    lines = path.read_text().splitlines()

    table_lines = [
        line for line in lines
        if line.strip().startswith("|")
    ]

    if len(table_lines) < 3:
        raise ValueError("No valid org table found.")

    # Header
    header = [
        col.strip()
        for col in table_lines[0].strip().strip("|").split("|")
    ]

    # Data rows
    data = []
    for line in table_lines[1:]:
        row = [
            col.strip()
            for col in line.strip().strip("|").split("|")
        ]

        if len(row) != len(header):
            continue

        data.append(row)

    df = pd.DataFrame(data, columns=header)

    # Convert numeric columns
    for col in df.columns:
        df[col] = pd.to_numeric(df[col])

    return df


def infer_h(n_cells: np.ndarray) -> np.ndarray:
    """
    Assume:
        h ~ 1 / sqrt(n_cells)

    which corresponds to uniform refinement
    in 2D square domain.
    """
    return 1.0 / np.sqrt(n_cells)


def fit_convergence_rate(h: np.ndarray, error: np.ndarray):
    """
    Fit:
        error = C * h^p

    using linear regression in log-log space.
    """

    log_h = np.log(h)
    log_e = np.log(error)

    slope, intercept = np.polyfit(log_h, log_e, 1)

    fitted = np.exp(intercept) * h ** slope

    return slope, fitted


def make_plot(
    h,
    error,
    quantity_name,
    source_path,
):
    slope, fitted = fit_convergence_rate(h, error)

    plt.figure(figsize=(6, 5))

    plt.loglog(
        h,
        error,
        "o",
        label=f"{quantity_name} data",
    )

    plt.loglog(
        h,
        fitted,
        "--",
        label=f"fit slope = {slope:.3f}",
    )

    plt.xlabel("h")
    plt.ylabel(quantity_name)
    plt.title(f"{quantity_name} convergence")

    plt.grid(True, which="both")
    plt.legend()

    output_path = (
        source_path.parent
        / f"{source_path.stem}_{quantity_name}.png"
    )

    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_path}")
    print(f"  slope = {slope:.6f}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot MMS convergence rates from dealii org table."
    )

    parser.add_argument(
        "table",
        type=Path,
        help="Path to org table file",
    )

    args = parser.parse_args()

    df = parse_org_table(args.table)

    h = infer_h(df["cells"].to_numpy())

    error_columns = [
        col for col in df.columns
        if col not in ["cycle", "cells", "dofs"]
    ]

    for col in error_columns:
        make_plot(
            h=h,
            error=df[col].values,
            quantity_name=col,
            source_path=args.table,
        )


if __name__ == "__main__":
    main()
