#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Funciones comunes para reproducir las figuras y tablas del manuscrito RMF."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

M = 1.0
OUTPUT_DIR = Path("output")
FIG_DIR = OUTPUT_DIR / "figures"
TABLE_DIR = OUTPUT_DIR / "tables"


def ensure_output_dirs():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def f_metric(r, mass=M):
    return 1.0 - 2.0 * mass / r


def v_eff(r, L2, kappa, mass=M):
    return f_metric(r, mass) * (kappa + L2 / r**2)


def dv_eff_dr(r, L2, kappa, mass=M):
    f = f_metric(r, mass)
    fp = 2.0 * mass / r**2
    return fp * (kappa + L2 / r**2) + f * (-2.0 * L2 / r**3)


def d2v_eff_dr2(r, L2, kappa, mass=M):
    return -4.0 * mass * kappa / r**3 + 6.0 * L2 / r**4 - 24.0 * mass * L2 / r**5


def circular_radii_timelike(L2, mass=M):
    disc = 1.0 - 12.0 * mass**2 / L2
    if disc < 0:
        return []
    if abs(disc) < 1e-14:
        return [L2 / (2.0 * mass)]
    root = math.sqrt(disc)
    return [
        (L2 / (2.0 * mass)) * (1.0 - root),
        (L2 / (2.0 * mass)) * (1.0 + root),
    ]


def b_crit(mass=M):
    return 3.0 * math.sqrt(3.0) * mass


def bisect_root(func, a, b, tol=1e-12, max_iter=500):
    fa, fb = func(a), func(b)
    if fa * fb > 0:
        raise ValueError("El intervalo no encierra una raíz.")
    for _ in range(max_iter):
        c = 0.5 * (a + b)
        fc = func(c)
        if 0.5 * abs(b - a) < tol:
            return c
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


@dataclass
class OrbitResult:
    x: np.ndarray
    y: np.ndarray
    r: np.ndarray
    phi: np.ndarray
    ur: np.ndarray
    max_constraint_error: float
    steps: int
    status: str


def rhs(state, E, L, kappa):
    _, r, ur, _ = state
    rr = max(r, 2.0 * M + 1e-12)
    return np.array([
        E / f_metric(rr),
        ur,
        -0.5 * dv_eff_dr(rr, L**2, kappa),
        L / rr**2,
    ], dtype=float)


def rk4_step(state, E, L, kappa, h):
    k1 = rhs(state, E, L, kappa)
    k2 = rhs(state + 0.5 * h * k1, E, L, kappa)
    k3 = rhs(state + 0.5 * h * k2, E, L, kappa)
    k4 = rhs(state + h * k3, E, L, kappa)
    return state + h * (k1 + 2*k2 + 2*k3 + k4) / 6.0


def radial_constraint_error(r, ur, E, L, kappa):
    return abs(ur**2 + v_eff(r, L**2, kappa) - E**2)


def integrate_orbit(E, L, kappa, r0, ur0, h=0.01, phi_goal=12*math.pi,
                    max_steps=1_000_000, r_stop=2.02, r_escape=None):
    state = np.array([0.0, r0, ur0, 0.0], dtype=float)
    xs, ys, rs, phis, urs = [], [], [], [], []
    max_err = 0.0
    status = "max_steps"

    for step in range(max_steps):
        _, r, ur, phi = state
        if not np.all(np.isfinite(state)):
            status = "non_finite"
            break
        if r <= r_stop:
            status = "horizon"
            break
        if phi >= phi_goal:
            status = "phi_goal"
            break
        if r_escape is not None and r >= r_escape and ur > 0 and step > 20:
            status = "escape"
            break

        xs.append(r * math.cos(phi))
        ys.append(r * math.sin(phi))
        rs.append(r)
        phis.append(phi)
        urs.append(ur)
        max_err = max(max_err, radial_constraint_error(r, ur, E, L, kappa))
        state = rk4_step(state, E, L, kappa, h)

    return OrbitResult(np.asarray(xs), np.asarray(ys), np.asarray(rs),
                       np.asarray(phis), np.asarray(urs), max_err,
                       len(rs), status)


def save_figure(fig, filename):
    ensure_output_dirs()
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, dpi=300, bbox_inches="tight")
