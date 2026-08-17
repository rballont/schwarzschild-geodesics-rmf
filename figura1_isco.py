#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figura 1: potencial efectivo e ISCO para L^2=12 M^2."""

import math
import numpy as np
import matplotlib.pyplot as plt

from schwarzschild_utils import M, integrate_orbit, save_figure, v_eff


def main():
    L2 = 12.0
    L = math.sqrt(L2)
    r0 = 6.0
    E = math.sqrt(v_eff(r0, L2, 1))

    orbit = integrate_orbit(E, L, 1, r0, 0.0, h=0.01, phi_goal=20*np.pi)
    r = np.linspace(2.05, 20.0, 2500)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(r/M, v_eff(r, L2, 1), lw=2, label=r"$V_{\rm eff}$")
    ax1.axhline(E**2, ls="--", lw=1.4, label=rf"$E^2={E**2:.4f}$")
    ax1.axvline(2, ls="--", lw=1.4, label=r"Horizonte $r=2M$")
    ax1.axvline(6, ls="--", lw=1.8, label=r"ISCO $r=6M$")
    ax1.scatter([6], [E**2], s=55, zorder=5)
    ax1.set_xlim(1.8, 20)
    ax1.set_xlabel(r"$r/M$")
    ax1.set_ylabel(r"$V_{\rm eff}, E^2$")
    ax1.set_title("Potencial efectivo")
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(fontsize=8)

    ax2.plot(orbit.x/M, orbit.y/M, lw=2, label="RK4")
    horizon = plt.Circle((0, 0), 2, fill=False, linestyle="--", linewidth=1.4)
    ax2.add_patch(horizon)
    ax2.set_aspect("equal", "box")
    ax2.set_xlabel(r"$x/M$")
    ax2.set_ylabel(r"$y/M$")
    ax2.set_title("Órbita circular ISCO")
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(fontsize=8)

    save_figure(fig, "figura1.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
