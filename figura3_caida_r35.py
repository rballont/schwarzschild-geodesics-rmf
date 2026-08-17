#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figura 3: caída directa para L^2=16 M^2 y r0=3.5 M."""

import math
import numpy as np
import matplotlib.pyplot as plt

from schwarzschild_utils import M, circular_radii_timelike, integrate_orbit, save_figure, v_eff


def circle(ax, radius, label):
    th = np.linspace(0, 2*np.pi, 500)
    ax.plot(radius*np.cos(th), radius*np.sin(th), lw=1.3, label=label)


def main():
    L2 = 16.0
    L = math.sqrt(L2)
    r0 = 3.5
    E = math.sqrt(v_eff(r0, L2, 1))
    r_minus, r_plus = circular_radii_timelike(L2)
    orbit = integrate_orbit(E, L, 1, r0, 0.0, h=0.002, phi_goal=6*np.pi,
                            r_stop=2.02, max_steps=500_000)
    r = np.linspace(2.05, 40.0, 3000)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(r/M, v_eff(r, L2, 1), lw=2, label=r"$V_{\rm eff}$")
    ax1.axhline(E**2, ls="--", lw=1.3, label=rf"$E^2(r_0)={E**2:.3f}$")
    ax1.axvline(2, ls=":", lw=1.6, label=r"Horizonte $r_h=2M$")
    ax1.axvline(r_minus, ls="--", lw=1.5, label=rf"$r_-={r_minus:.1f}M$")
    ax1.axvline(r_plus, ls="--", lw=1.5, label=rf"$r_+={r_plus:.1f}M$")
    ax1.axvline(r0, ls=":", lw=1.5, label=rf"$r_0={r0:.1f}M$")
    ax1.axvspan(2, r_minus, alpha=0.12, label="Región I")
    ax1.axvspan(r_minus, r_plus, alpha=0.16, label="Región II")
    ax1.set_xlim(1.8, 40)
    ax1.set_xlabel(r"$r/M$")
    ax1.set_ylabel(r"$V_{\rm eff}, E^2$")
    ax1.set_title(r"Potencial efectivo, $L^2=16M^2$")
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(fontsize=7, ncol=2)

    ax2.plot(orbit.x/M, orbit.y/M, lw=2, label=rf"Trayectoria RK4, $r_0={r0:.1f}M$")
    horizon = plt.Circle((0, 0), 2, fill=False, linestyle="--", linewidth=1.4)
    ax2.add_patch(horizon)
    circle(ax2, r_minus, rf"$r_-={r_minus:.1f}M$")
    circle(ax2, r_plus, rf"$r_+={r_plus:.1f}M$")
    ax2.set_aspect("equal", "box")
    ax2.set_xlabel(r"$x/M$")
    ax2.set_ylabel(r"$y/M$")
    ax2.set_title("Trayectoria de caída directa")
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.32))

    save_figure(fig, "figura3.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
