#!/usr/bin/env python3
"""
Visualize the manufactured solution and its source terms on [0, 1] x [0, 1].

Produces a 2x3 figure:
  - velocity vector field
  - pressure
  - temperature
  - momentum RHS x-component
  - momentum RHS y-component
  - temperature RHS
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


def build_functions():
    x, y, nu, kappa = sp.symbols("x y nu kappa", real=True)
    pi = sp.pi

    u_x = sp.sin(pi * x) * sp.cos(pi * y)
    u_y = -sp.cos(pi * x) * sp.sin(pi * y)
    p = sp.cos(pi * x) * sp.sin(pi * y)
    T = sp.cos(2 * pi * x) * sp.sin(pi * y)

    grad_ux = (sp.diff(u_x, x), sp.diff(u_x, y))
    grad_uy = (sp.diff(u_y, x), sp.diff(u_y, y))
    lap_ux = sp.diff(u_x, x, 2) + sp.diff(u_x, y, 2)
    lap_uy = sp.diff(u_y, x, 2) + sp.diff(u_y, y, 2)
    momentum_rhs_x = sp.simplify(u_x * grad_ux[0] + u_y * grad_ux[1] - nu * lap_ux + sp.diff(p, x))
    momentum_rhs_y = sp.simplify(u_x * grad_uy[0] + u_y * grad_uy[1] - nu * lap_uy + sp.diff(p, y))

    grad_T = (sp.diff(T, x), sp.diff(T, y))
    lap_T = sp.diff(T, x, 2) + sp.diff(T, y, 2)
    temperature_rhs = sp.simplify(u_x * grad_T[0] + u_y * grad_T[1] - kappa * lap_T)

    funcs = {
        "u_x": sp.lambdify((x, y), u_x, "numpy"),
        "u_y": sp.lambdify((x, y), u_y, "numpy"),
        "p": sp.lambdify((x, y), p, "numpy"),
        "T": sp.lambdify((x, y), T, "numpy"),
        "fx": sp.lambdify((x, y, nu), momentum_rhs_x, "numpy"),
        "fy": sp.lambdify((x, y, nu), momentum_rhs_y, "numpy"),
        "s": sp.lambdify((x, y, kappa), temperature_rhs, "numpy"),
    }
    return funcs


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nu", type=float, default=1.0, help="Kinematic viscosity used in the momentum RHS.")
    parser.add_argument("--kappa", type=float, default=1.0, help="Thermal diffusivity used in the temperature RHS.")
    parser.add_argument("--n", type=int, default=200, help="Grid resolution per axis.")
    parser.add_argument("--output", type=Path, default=Path("mms_visualization.png"))
    parser.add_argument("--show", action="store_true", help="Show the figure interactively.")
    return parser.parse_args()


def main():
    args = parse_args()
    f = build_functions()

    x = np.linspace(0.0, 1.0, args.n)
    y = np.linspace(0.0, 1.0, args.n)
    X, Y = np.meshgrid(x, y, indexing="xy")

    U = f["u_x"](X, Y)
    V = f["u_y"](X, Y)
    P = f["p"](X, Y)
    T = f["T"](X, Y)
    FX = f["fx"](X, Y, args.nu)
    FY = f["fy"](X, Y, args.nu)
    S = f["s"](X, Y, args.kappa)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    fig.suptitle("Manufactured solution and source terms", fontsize=16)

    def contour(ax, Z, title, cmap="coolwarm"):
        im = ax.contourf(X, Y, Z, levels=40, cmap=cmap)
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, shrink=0.9)

    ax = axes[0, 0]
    ax.quiver(X[::10, ::10], Y[::10, ::10], U[::10, ::10], V[::10, ::10], color="black", pivot="mid", scale=25)
    ax.set_title("Velocity field")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    contour(axes[0, 1], P, "Pressure")
    contour(axes[0, 2], T, "Temperature", cmap="viridis")
    contour(axes[1, 0], FX, r"Momentum RHS $f_x$")
    contour(axes[1, 1], FY, r"Momentum RHS $f_y$")
    contour(axes[1, 2], S, r"Temperature RHS $s$", cmap="magma")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    print(f"saved {args.output}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
