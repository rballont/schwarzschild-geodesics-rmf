#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figura 2: potencial efectivo para L^2=16 M^2."""

import numpy as np
import matplotlib.pyplot as plt

from schwarzschild_utils import M, circular_radii_timelike, save_figure, v_eff


def main():
    L2 = 16.0
    r_minus, r_plus = circular_radii_timelike(L2)
    r = np.linspace(2.05, 30.0, 3000)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(r/M, v_eff(r, L2, 1), lw=2, label=r"$V_{\rm eff}(r)$")
    ax.axhline(1, ls=":", lw=1.3, label=r"$E^2=1$")
    ax.axvspan(2, r_minus, alpha=0.12, label="Región I: caída directa")
    ax.axvspan(r_minus, r_plus, alpha=0.16, label="Región II: órbitas posibles")
    ax.axvline(2, ls=":", lw=1.6, label=r"Horizonte $r_h=2M$")
    ax.axvline(r_minus, ls="--", lw=1.6, label=rf"$r_-={r_minus:.1f}M$")
    ax.axvline(r_plus, ls="--", lw=1.6, label=rf"$r_+={r_plus:.1f}M$")
    ax.scatter([r_minus, r_plus], [v_eff(r_minus, L2, 1), v_eff(r_plus, L2, 1)], s=55, zorder=5)
    ax.set_xlim(1.8, 30)
    ax.set_xlabel(r"$r/M$")
    ax.set_ylabel(r"$V_{\rm eff}(r)$")
    ax.set_title(r"Potencial efectivo tipo tiempo, $L^2=16M^2$")
    ax.grid(True, ls=":", alpha=0.35)
    ax.legend(fontsize=8, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.36))

    save_figure(fig, "figura2.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
