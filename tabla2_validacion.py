#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tabla II: comparación entre valores analíticos y raíces numéricas."""

import math

from schwarzschild_utils import M, TABLE_DIR, b_crit, bisect_root, d2v_eff_dr2, dv_eff_dr, ensure_output_dirs, f_metric


def analytical_radii(L2):
    disc = math.sqrt(1.0 - 12.0 / L2)
    return 0.5 * L2 * (1.0 - disc), 0.5 * L2 * (1.0 + disc)


def L2_circular(r):
    return M * r**2 / (r - 3.0 * M)


def main():
    ensure_output_dirs()

    r_isco = bisect_root(lambda r: d2v_eff_dr2(r, L2_circular(r), 1), 5.0, 7.0)
    r16m = bisect_root(lambda r: dv_eff_dr(r, 16.0, 1), 2.01, 5.99)
    r16p = bisect_root(lambda r: dv_eff_dr(r, 16.0, 1), 6.01, 40.0)
    r24m = bisect_root(lambda r: dv_eff_dr(r, 24.0, 1), 2.01, 5.99)
    r24p = bisect_root(lambda r: dv_eff_dr(r, 24.0, 1), 6.01, 40.0)
    rph = bisect_root(lambda r: dv_eff_dr(r, 36.0, 0), 2.01, 10.0)
    bc_num = rph / math.sqrt(f_metric(rph))

    a16m, a16p = analytical_radii(16.0)
    a24m, a24p = analytical_radii(24.0)

    rows = [
        ("ISCO", 6.0, r_isco),
        ("r-, L2=16M2", a16m, r16m),
        ("r+, L2=16M2", a16p, r16p),
        ("r-, L2=24M2", a24m, r24m),
        ("r+, L2=24M2", a24p, r24p),
        ("r_ph", 3.0, rph),
        ("b_crit", b_crit(), bc_num),
    ]

    lines = ["Magnitud | Analítico | Numérico | Error absoluto\n"]
    for name, analytic, numeric in rows:
        lines.append(f"{name} | {analytic:.10f} | {numeric:.10f} | {abs(numeric-analytic):.3e}\n")

    out = TABLE_DIR / "tabla2_validacion.txt"
    out.write_text("".join(lines), encoding="utf-8")
    print("".join(lines))


if __name__ == "__main__":
    main()
