#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figura 5: potencial efectivo para geodésicas nulas."""

import numpy as np
import matplotlib.pyplot as plt

from schwarzschild_utils import M, save_figure, v_eff


def main():
    L2 = 36.0
    r = np.linspace(2.05, 25.0, 2500)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(r/M, v_eff(r, L2, 0), lw=2, label=r"$V_{\rm eff}$")
    ax.axvline(2, ls="--", lw=1.5, label=r"$r_h=2M$")
    ax.axvline(3, ls="--", lw=1.7, label=r"$r=3M$")
    ax.scatter([3], [v_eff(3, L2, 0)], s=55, zorder=5)
    ax.set_xlabel(r"$r/M$")
    ax.set_ylabel(r"$V_{\rm eff}$")
    ax.set_title(r"Potencial efectivo nulo, $L^2=36M^2$")
    ax.grid(True, ls=":", alpha=0.35)
    ax.legend(fontsize=8)

    save_figure(fig, "figura5.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
