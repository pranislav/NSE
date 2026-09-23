#!/usr/bin/env python3

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


MAX_TRAILING_OUTLIERS = 2
MIN_FIT_POINTS = 3
OUTLIER_MIN_LOG_DEVIATION = math.log(1.25)
OUTLIER_MIN_RATE_DROP = 0.75
OUTLIER_RATE_STD_MULTIPLIER = 3


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


def fit_convergence_rate(
    h: np.ndarray,
    error: np.ndarray,
    evaluation_h: np.ndarray | None = None,
):
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

    if evaluation_h is None:
        evaluation_h = h
    fitted = np.exp(intercept) * evaluation_h ** slope

    slope_std = np.sqrt(cov[0, 0])

    return slope, fitted, slope_std


def find_trailing_outliers(h: np.ndarray, error: np.ndarray) -> np.ndarray:
    """Identify up to two fine-grid errors that depart from the coarse-grid trend.

    Only the final refinement levels are considered because round-off effects are
    expected at the smallest errors. A point is excluded only when it is at least
    25% above the extrapolated error *and* its local rate drops substantially
    relative to the preceding refinement levels. This prevents smooth changes in
    slope from being treated as outliers.
    """
    n_points = len(error)
    max_candidates = min(MAX_TRAILING_OUTLIERS, n_points - MIN_FIT_POINTS)

    for n_outliers in range(max_candidates, 0, -1):
        retained = slice(0, n_points - n_outliers)
        log_h = np.log(h[retained])
        log_error = np.log(error[retained])
        slope, intercept = np.polyfit(log_h, log_error, deg=1)

        excluded = np.arange(n_points - n_outliers, n_points)
        extrapolated_log_error = slope * np.log(h[excluded]) + intercept
        deviations = np.log(error[excluded]) - extrapolated_log_error

        local_rates = np.log(error[:-1] / error[1:]) / np.log(h[:-1] / h[1:])
        candidate_local_rates = local_rates[n_points - n_outliers - 1:]
        reference_local_rates = local_rates[:n_points - n_outliers - 1]
        reference_rate = np.median(reference_local_rates)
        reference_spread = np.std(reference_local_rates, ddof=1)
        rate_drop_threshold = max(
            OUTLIER_MIN_RATE_DROP,
            OUTLIER_RATE_STD_MULTIPLIER * reference_spread,
        )
        rate_drops = reference_rate - candidate_local_rates

        if (
            np.all(deviations > OUTLIER_MIN_LOG_DEVIATION)
            and np.all(rate_drops > rate_drop_threshold)
        ):
            return excluded

    return np.array([], dtype=int)


def make_plot_and_log_rates(
    h,
    error,
    quantity_name,
    source_path,
    rates_file,
    outlier_indices,
):
    retained = np.ones(len(error), dtype=bool)
    retained[outlier_indices] = False
    slope, fitted, slope_std = fit_convergence_rate(
        h[retained],
        error[retained],
        evaluation_h=h,
    )

    plt.figure(figsize=(6, 5))

    plt.loglog(
        h[retained],
        error[retained],
        "o",
        label=f"{quantity_name} data",
    )

    if len(outlier_indices):
        plt.loglog(
            h[outlier_indices],
            error[outlier_indices],
            "o",
            color="red",
            label="excluded outlier",
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
    if len(outlier_indices):
        cycles = ", ".join(str(index) for index in outlier_indices)
        rates_file.write(f"  excluded outlier cycle(s): {cycles}\n")
    print(f"Estimated convergence rate for {quantity_name}: {rate_with_uncertainty}")
    if len(outlier_indices):
        print(f"Excluded {quantity_name} outlier cycle(s): {cycles}")
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


def format_outlier_tex(value: str, is_outlier: bool) -> str:
    if is_outlier:
        return rf"\textcolor{{red}}{{{value}}}"
    return value


def infer_pressure_degree(source_path: Path) -> int:
    """Infer the pressure polynomial degree from the input table name or path."""
    for text in (source_path.stem, *source_path.parts[::-1]):
        match = re.search(r"(?:^|[_-])q(?P<degree>\d+)(?:$|[_-])", text)
        if match:
            return int(match.group("degree"))

        match = re.search(r"(?:^|[_-])deg(?P<degree>\d+)(?:$|[_-])", text)
        if match:
            return int(match.group("degree"))

    raise ValueError(
        "Could not infer the polynomial degree from the table path. "
        "Expected a name containing q<N> or deg<N>."
    )


def theoretical_convergence_rate(error_column: str, pressure_degree: int) -> int:
    """Return the expected rate, accounting for the higher velocity/temperature degree."""
    basis_degree = pressure_degree
    if error_column.endswith(("_velocity", "_temperature")):
        basis_degree += 1

    if error_column.startswith("H1_"):
        return basis_degree
    if error_column.startswith("L2_"):
        return basis_degree + 1

    raise ValueError(
        f"Cannot determine the theoretical convergence rate for {error_column!r}."
    )


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


def write_local_convergence_table(
    df,
    h,
    error_columns,
    fitted_rates,
    theoretical_rates,
    pressure_degree,
    outliers,
    output_path,
):
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
        r"% Requires \usepackage{xcolor}",
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
            is_outlier = row_index in outliers[col]
            cells.append(format_outlier_tex(format_error(row[col]), is_outlier))
            if row_index == 0:
                cells.append("--")
            else:
                cells.append(format_outlier_tex(
                    format_local_rate(local_rates[col][row_index]), is_outlier
                ))

        lines.append(" & ".join(cells) + r" \\ \hline")

    fitted_row = [
        rf"\multicolumn{{3}}{{|r|}}{{fitted rate}}",
    ]
    for col in error_columns:
        fitted_row.extend(["--", format_tex_rate_with_uncertainty(fitted_rates[col])])
    lines.append(" & ".join(fitted_row) + r" \\ \hline")

    theoretical_row = [
        rf"\multicolumn{{3}}{{|r|}}{{theoretical rate}}",
    ]
    for col in error_columns:
        theoretical_row.extend(["--", f"${theoretical_rates[col]}$"])
    lines.append(" & ".join(theoretical_row) + r" \\ \hline")

    lines.extend([
        r"\end{tabular}",
        rf"\caption{{Local convergence rates computed from adjacent refinement levels. Pressure is approximated with degree-${pressure_degree}$ polynomials, while velocity and temperature are approximated with degree-${pressure_degree + 1}$ polynomials. Red values are fine-grid outliers excluded from the log-log least-squares fit. The final two rows show fitted and theoretical rates, respectively.}}",
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

    pressure_degree = infer_pressure_degree(args.table)
    theoretical_rates = {
        col: theoretical_convergence_rate(col, pressure_degree)
        for col in error_columns
    }
    outliers = {
        col: find_trailing_outliers(h, df[col].to_numpy())
        for col in error_columns
    }

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
                outlier_indices=outliers[col],
            )

    write_local_convergence_table(
        df=df,
        h=h,
        error_columns=error_columns,
        fitted_rates=fitted_rates,
        theoretical_rates=theoretical_rates,
        pressure_degree=pressure_degree,
        outliers=outliers,
        output_path=args.table.parent / "convergence_rates.tex",
    )


if __name__ == "__main__":
    main()
