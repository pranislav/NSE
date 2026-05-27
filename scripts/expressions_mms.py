"""
Shared manufactured-solution expressions for the 2D conjugate heat transfer
project.

This is the single source of truth for the exact fields used by the other
scripts:
  - velocity u = (sin(pi x) cos(pi y), -cos(pi x) sin(pi y))
  - pressure p = cos(pi x) sin(pi y)
  - temperature T = cos(2 pi x) sin(pi y)

The module also builds the corresponding source terms for:
  - steady incompressible Navier-Stokes
  - steady advection-diffusion temperature equation
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ExprSet:
    u: tuple[sp.Expr, sp.Expr]
    p: sp.Expr
    T: sp.Expr
    grad_u: list[list[sp.Expr]]
    grad_p: tuple[sp.Expr, sp.Expr]
    grad_T: tuple[sp.Expr, sp.Expr]
    momentum_rhs: tuple[sp.Expr, sp.Expr]
    continuity_rhs: sp.Expr
    temperature_rhs: sp.Expr


def build_expressions():
    x, y, nu, thermal_diffusivity = sp.symbols("x y nu thermal_diffusivity", real=True)
    pi = sp.pi

    u = (
        sp.sin(pi * x) * sp.cos(pi * y),
        -sp.cos(pi * x) * sp.sin(pi * y),
    )
    p = sp.cos(pi * x) * sp.sin(pi * y)
    T = sp.cos(2 * pi * x) * sp.sin(pi * y)

    ux, uy = u
    grad_u = [[sp.diff(ux, x), sp.diff(ux, y)], [sp.diff(uy, x), sp.diff(uy, y)]]
    lap_u = [sp.diff(ux, x, 2) + sp.diff(ux, y, 2), sp.diff(uy, x, 2) + sp.diff(uy, y, 2)]
    grad_p = (sp.diff(p, x), sp.diff(p, y))
    adv_u = (ux * grad_u[0][0] + uy * grad_u[0][1], ux * grad_u[1][0] + uy * grad_u[1][1])
    momentum_rhs = tuple(sp.simplify(adv_u[i] - nu * lap_u[i] + grad_p[i]) for i in range(2))
    continuity_rhs = sp.simplify(sp.diff(ux, x) + sp.diff(uy, y))

    grad_T = (sp.diff(T, x), sp.diff(T, y))
    lap_T = sp.diff(T, x, 2) + sp.diff(T, y, 2)
    adv_T = ux * grad_T[0] + uy * grad_T[1]
    temperature_rhs = sp.simplify(adv_T - thermal_diffusivity * lap_T)

    return ExprSet(
        u=u,
        p=p,
        T=T,
        grad_u=grad_u,
        grad_p=grad_p,
        grad_T=grad_T,
        momentum_rhs=momentum_rhs,
        continuity_rhs=continuity_rhs,
        temperature_rhs=temperature_rhs,
    ), (x, y, nu, thermal_diffusivity)


def boundary_values(expr: sp.Expr, x: sp.Symbol, y: sp.Symbol):
    return {
        "x = 0": sp.simplify(expr.subs(x, 0)),
        "x = 1": sp.simplify(expr.subs(x, 1)),
        "y = 0": sp.simplify(expr.subs(y, 0)),
        "y = 1": sp.simplify(expr.subs(y, 1)),
    }
