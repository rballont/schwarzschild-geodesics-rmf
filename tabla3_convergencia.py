#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tabla III: conservación de la restricción radial para distintos pasos RK4."""

import math

from schwarzschild_utils import TABLE_DIR, ensure_output_dirs, integrate_orbit, v_eff


def main():
    ensure_output_dirs()

    E = 1.0
    b = 7.0
    L = b * E
    r0 = 20.0
    ur0 = -math.sqrt(E**2 - v_eff(r0, L**2, 0))

    rows = []
    for h in [0.04, 0.02, 0.01, 0.005]:
        res = integrate_orbit(E, L, 0, r0, ur0, h=h,
                              phi_goal=4*math.pi,
                              r_escape=30.0,
                              max_steps=500_000)
        rows.append((h, res.max_constraint_error, res.steps, res.status))

    lines = ["h | Delta_max | pasos | estado\n"]
    for h, err, steps, status in rows:
        lines.append(f"{h:.3f} | {err:.3e} | {steps} | {status}\n")

    out = TABLE_DIR / "tabla3_convergencia.txt"
    out.write_text("".join(lines), encoding="utf-8")
    print("".join(lines))


if __name__ == "__main__":
    main()
