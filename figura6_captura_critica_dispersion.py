#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figura 6: captura, trayectoria crítica y dispersión de geodésicas nulas."""

import math
import numpy as np
import matplotlib.pyplot as plt

from schwarzschild_utils import M, b_crit, f_metric, rk4_step, v_eff, save_figure

YELLOW = "#F2C94C"
BLACK = "#000000"


def integrate_null_case(factor, h=0.002, r0=30.0, max_steps=3_000_000):
    E = 1.0
    L = factor * b_crit()
    ur0 = -math.sqrt(E**2 - v_eff(r0, L**2, 0))
    state = np.array([0.0, r0, ur0, 0.0], dtype=float)

    xs, ys = [], []
    for step in range(max_steps):
        _, r, ur, phi = state
        if not np.all(np.isfinite(state)):
            break

        xs.append(r * math.cos(phi))
        ys.append(r * math.sin(phi))

        if factor < 1.0 and r <= 2.02:
            break
        if abs(factor - 1.0) < 1e-14 and r <= 3.00005 and abs(ur) < 5e-4:
            break
        if factor > 1.0 and r >= r0 and ur > 0 and step > 100:
            break
        if phi > 20.0 * math.pi:
            break

        state = rk4_step(state, E, L, 0, h)

    return np.asarray(xs), np.asarray(ys)


def main():
    cases = [
        (0.9, r"(a) Captura: $b=0.9\,b_{\rm crit}$", (-6, 31), (-8, 12)),
        (1.0, r"(b) Trayectoria crítica: $b=b_{\rm crit}$", (-10, 31), (-12, 12)),
        (1.2, r"(c) Dispersión: $b=1.2\,b_{\rm crit}$", (-8, 31), (-18, 14)),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    for ax, (factor, title, xlim, ylim) in zip(axes, cases):
        x, y = integrate_null_case(factor)
        ax.plot(x/M, y/M, lw=2.1, color=YELLOW, label="Geodésica nula")

        horizon = plt.Circle((0, 0), 2.0, facecolor=BLACK, edgecolor=BLACK,
                             label=r"Horizonte $r=2M$")
        photon_sphere = plt.Circle((0, 0), 3.0, fill=False, edgecolor=YELLOW,
                                   linestyle="--", linewidth=1.5,
                                   label=r"Esfera de fotones $r=3M$")
        ax.add_patch(horizon)
        ax.add_patch(photon_sphere)

        ax.set_aspect("equal", "box")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"$x/M$")
        ax.set_ylabel(r"$y/M$")
        ax.set_title(title)
        ax.grid(True, ls=":", alpha=0.30)
        ax.legend(fontsize=8, loc="best")

    save_figure(fig, "figura6.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
