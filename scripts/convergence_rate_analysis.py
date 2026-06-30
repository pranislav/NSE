#!/usr/bin/env python3

import argparse
import math
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

    coef, cov = np.polyfit(log_h, log_e, deg=1, cov=True)

    slope = coef[0]
    intercept = coef[1]

    fitted = np.exp(intercept) * h ** slope

    slope_std = np.sqrt(cov[0, 0])

    return slope, fitted, slope_std


def make_plot_and_log_rates(
    h,
    error,
    quantity_name,
    source_path,
    rates_file
):
    slope, fitted, slope_std = fit_convergence_rate(h, error)

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
    rate_with_uncertainty = format_with_uncertainty(slope, slope_std)
    rates_file.write(f"{quantity_name}: {rate_with_uncertainty}\n")
    print(f"Estimated convergence rate for {quantity_name}: {rate_with_uncertainty}")
    return rate_with_uncertainty


def tex_escape(text):
    return (
        str(text)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def format_error(value):
    return f"{value:.3e}"


def format_local_rate(value):
    if not np.isfinite(value):
        return "--"
    return f"{value:.2f}"


def format_tex_rate_with_uncertainty(rate):
    tex_rate = tex_escape(rate).replace("±", r"\pm")
    return f"${tex_rate}$"


def compute_local_rates(h, error):
    rates = [None]
    for previous_h, current_h, previous_error, current_error in zip(
        h[:-1],
        h[1:],
        error[:-1],
        error[1:],
    ):
        if previous_h <= 0 or current_h <= 0 or previous_error <= 0 or current_error <= 0:
            rates.append(np.nan)
            continue

        rates.append(
            np.log(previous_error / current_error)
            / np.log(previous_h / current_h)
        )

    return rates


def write_local_convergence_table(df, h, error_columns, fitted_rates, output_path):
    local_rates = {
        col: compute_local_rates(h, df[col].to_numpy())
        for col in error_columns
    }

    column_format = "|r|r|r|" + "r|r|" * len(error_columns)
    headers = ["cycle", r"\# cells", r"\# dofs"]
    for col in error_columns:
        quantity = tex_escape(col)
        headers.extend([quantity, "rate"])

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\begin{{tabular}}{{{column_format}}}",
        r"\hline",
        " & ".join(headers) + r" \\ \hline",
    ]

    for row_index, row in df.iterrows():
        cells = [
            f"{int(row['cycle'])}",
            f"{int(row['cells'])}",
            f"{int(row['dofs'])}",
        ]

        for col in error_columns:
            cells.append(format_error(row[col]))
            if row_index == 0:
                cells.append("--")
            else:
                cells.append(format_local_rate(local_rates[col][row_index]))

        lines.append(" & ".join(cells) + r" \\ \hline")

    fitted_row = [
        rf"\multicolumn{{3}}{{|r|}}{{fitted rate}}",
    ]
    for col in error_columns:
        fitted_row.extend(["--", format_tex_rate_with_uncertainty(fitted_rates[col])])
    lines.append(" & ".join(fitted_row) + r" \\ \hline")

    lines.extend([
        r"\end{tabular}",
        r"\caption{Local convergence rates computed from adjacent refinement levels. The final row shows fitted rates with uncertainties from the log-log least-squares fit.}",
        r"\end{table}",
        "",
    ])

    output_path.write_text("\n".join(lines))
    print(f"Saved: {output_path}")

def format_with_uncertainty(value, error):
    if error <= 0:
        return f"{value}"

    exponent = math.floor(math.log10(abs(error)))
    first_digit = error / 10**exponent

    # 2 significant digits if leading digit is 1 or 2
    sig_digits = 2 if first_digit < 3 else 1

    decimals = -(exponent - (sig_digits - 1))

    value_rounded = round(value, decimals)
    error_rounded = round(error, decimals)

    return f"{value_rounded:.{max(0, decimals)}f} ± {error_rounded:.{max(0, decimals)}f}"

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

    fitted_rates = {}
    rates_path = args.table.parent / "convergence_rates.txt"
    with rates_path.open("w") as rates_file:
        for col in error_columns:
            fitted_rates[col] = make_plot_and_log_rates(
                h=h,
                error=df[col].values,
                quantity_name=col,
                source_path=args.table,
                rates_file=rates_file,
            )

    write_local_convergence_table(
        df=df,
        h=h,
        error_columns=error_columns,
        fitted_rates=fitted_rates,
        output_path=args.table.parent / "convergence_rates.tex",
    )


if __name__ == "__main__":
    main()
