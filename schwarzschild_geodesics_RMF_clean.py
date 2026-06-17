#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
schwarzschild_geodesics_RMF_clean.py

"Clasificación analítico-computacional de geodésicas tipo tiempo y nulas en el
espacio-tiempo de Schwarzschild".

El archivo reorganiza los códigos del Anexo D de la tesis y genera las figuras
principales del manuscrito para RMF, además de tablas de validación numérica.

Salidas:
    output_rmf/figures/figura1.png ... figura6.png
    output_rmf/tables/tabla_validacion_RMF.txt
    output_rmf/tables/tabla_validacion_RMF.tex

Requisitos:
    pip install numpy matplotlib

Autores:
    Ricardo Angelo Ballon Tito
    Rolando Moisés Perca Gonzales
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Configuración general
# ============================================================

M = 1.0
OUTPUT_DIR = Path("output_rmf")
FIG_DIR = OUTPUT_DIR / "figures"
TAB_DIR = OUTPUT_DIR / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Potenciales y radios característicos
# ============================================================

def f_metric(r, mass=M):
    """Función métrica de Schwarzschild: f(r)=1-2M/r."""
    return 1.0 - 2.0 * mass / r


def v_eff(r, L2, kappa, mass=M):
    """
    Potencial efectivo:
        V_eff(r) = (1 - 2M/r) (kappa + L^2/r^2).

    kappa = 1: geodésicas tipo tiempo.
    kappa = 0: geodésicas nulas.
    """
    return f_metric(r, mass) * (kappa + L2 / r**2)


def dv_eff_dr(r, L2, kappa, mass=M):
    """Derivada analítica de V_eff respecto a r."""
    f = f_metric(r, mass)
    fp = 2.0 * mass / r**2
    return fp * (kappa + L2 / r**2) + f * (-2.0 * L2 / r**3)


def circular_radii_timelike(L2, mass=M):
    """
    Radios críticos tipo tiempo:
        r_± = L^2/(2M) [1 ± sqrt(1 - 12M^2/L^2)].
    """
    disc = 1.0 - 12.0 * mass**2 / L2
    if disc < 0:
        return []
    if abs(disc) < 1e-14:
        return [L2 / (2.0 * mass)]
    root = math.sqrt(disc)
    r_minus = (L2 / (2.0 * mass)) * (1.0 - root)
    r_plus = (L2 / (2.0 * mass)) * (1.0 + root)
    return [r_minus, r_plus]


def photon_sphere_radius(mass=M):
    return 3.0 * mass


def b_crit(mass=M):
    return 3.0 * math.sqrt(3.0) * mass


def impact_parameter_from_periastron(rp, mass=M):
    """b(r_p)=r_p/sqrt(1-2M/r_p), válido para r_p>2M."""
    return rp / np.sqrt(f_metric(rp, mass))


# ============================================================
# Integración RK4
# ============================================================

@dataclass
class OrbitResult:
    x: np.ndarray
    y: np.ndarray
    r: np.ndarray
    phi: np.ndarray
    ur: np.ndarray
    max_constraint_error: float
    final_constraint_error: float
    steps: int
    status: str


def rhs(state, E, L, kappa):
    """
    Sistema de primer orden:
        dt/dlambda   = E/f(r)
        dr/dlambda   = u_r
        du_r/dlambda = -1/2 dV_eff/dr
        dphi/dlambda = L/r^2
    """
    _, r, ur, _ = state
    if r <= 2.0 * M:
        r = 2.0 * M + 1e-12
    return np.array([
        E / f_metric(r),
        ur,
        -0.5 * dv_eff_dr(r, L**2, kappa),
        L / r**2,
    ], dtype=float)


def rk4_step(state, E, L, kappa, h):
    k1 = rhs(state, E, L, kappa)
    k2 = rhs(state + 0.5 * h * k1, E, L, kappa)
    k3 = rhs(state + 0.5 * h * k2, E, L, kappa)
    k4 = rhs(state + h * k3, E, L, kappa)
    return state + (h / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)


def radial_constraint_error(r, ur, E, L, kappa):
    """Delta = |u_r^2 + V_eff(r) - E^2|."""
    if r <= 2.0 * M:
        return float("nan")
    return abs(ur**2 + v_eff(r, L**2, kappa) - E**2)


def integrate_orbit(E, L, kappa, r0, ur0, h=0.01, phi_goal=12*math.pi,
                    max_steps=1_000_000, r_stop=2.02, r_escape=None):
    """Integra una trayectoria en el plano ecuatorial."""
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

        err = radial_constraint_error(r, ur, E, L, kappa)
        if np.isfinite(err):
            max_err = max(max_err, err)
        state = rk4_step(state, E, L, kappa, h)

    _, r, ur, _ = state
    final_err = radial_constraint_error(r, ur, E, L, kappa)
    return OrbitResult(np.array(xs), np.array(ys), np.array(rs), np.array(phis),
                       np.array(urs), max_err, final_err, len(rs), status)


# ============================================================
# Utilidades gráficas
# ============================================================

def draw_horizon(ax, radius=2.0*M, label=r"Horizonte $r=2M$"):
    th = np.linspace(0, 2*np.pi, 500)
    ax.plot(radius*np.cos(th), radius*np.sin(th), "k--", lw=1.4, label=label)


def draw_circle(ax, radius, linestyle="-", label=None, lw=1.3):
    th = np.linspace(0, 2*np.pi, 500)
    ax.plot(radius*np.cos(th), radius*np.sin(th), linestyle, lw=lw, label=label)


def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(FIG_DIR / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Figuras RMF
# ============================================================

def figure1_isco():
    """Figura 1: ISCO para L^2=12M^2."""
    L2 = 12.0
    L = math.sqrt(L2)
    r0 = 6.0
    E = math.sqrt(v_eff(r0, L2, 1))
    orbit = integrate_orbit(E, L, 1, r0, 0.0, h=0.02, phi_goal=20*np.pi)
    r_plot = np.linspace(2.05, 20, 2500)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(r_plot/M, v_eff(r_plot, L2, 1), lw=2, label=r"$V_{\rm eff}$")
    ax1.axhline(E**2, ls="--", lw=1.4, label=rf"$E^2={E**2:.4f}$")
    ax1.axvline(2, ls="--", lw=1.4, label=r"Horizonte $r=2M$")
    ax1.axvline(6, ls="--", lw=1.8, label=r"ISCO $r=6M$")
    ax1.scatter([6], [E**2], s=55, zorder=5, label=r"$r_0=6M$")
    ax1.set_xlim(1.8, 20)
    ax1.set_xlabel(r"$r/M$")
    ax1.set_ylabel(r"$V_{\rm eff},\,E^2$")
    ax1.set_title("Potencial efectivo")
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(fontsize=8)

    ax2.plot(orbit.x/M, orbit.y/M, lw=2, label="RK4")
    draw_horizon(ax2)
    ax2.set_aspect("equal", "box")
    ax2.set_xlabel(r"$x/M$")
    ax2.set_ylabel(r"$y/M$")
    ax2.set_title("Órbita circular ISCO")
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(fontsize=8)
    savefig(fig, "figura1.png")


def figure2_l16_potential():
    """Figura 2: potencial efectivo para L^2=16M^2."""
    L2 = 16.0
    r_minus, r_plus = circular_radii_timelike(L2)
    r_plot = np.linspace(2.05, 30, 3000)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(r_plot/M, v_eff(r_plot, L2, 1), lw=2, label=r"$V_{\rm eff}(r)$")
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
    savefig(fig, "figura2.png")


def figure3_direct_infall():
    """Figura 3: caída directa para L^2=16M^2 y r0=3.5M."""
    L2 = 16.0
    L = math.sqrt(L2)
    r0 = 3.5
    E = math.sqrt(v_eff(r0, L2, 1))
    orbit = integrate_orbit(E, L, 1, r0, 0.0, h=0.002, phi_goal=6*np.pi,
                            r_stop=2.02, max_steps=500_000)
    r_minus, r_plus = circular_radii_timelike(L2)
    r_plot = np.linspace(2.05, 40, 3000)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(r_plot/M, v_eff(r_plot, L2, 1), lw=2, label=r"$V_{\rm eff}$")
    ax1.axhline(E**2, ls="--", lw=1.3, label=rf"$E^2(r_0)={E**2:.3f}$")
    ax1.axvline(2, ls=":", lw=1.6, label=r"Horizonte $r_h=2M$")
    ax1.axvline(r_minus, ls="--", lw=1.5, label=rf"$r_-={r_minus:.1f}M$")
    ax1.axvline(r_plus, ls="--", lw=1.5, label=rf"$r_+={r_plus:.1f}M$")
    ax1.axvline(r0, ls=":", lw=1.5, label=rf"$r_0={r0:.1f}M$")
    ax1.axvspan(2, r_minus, alpha=0.12, label="Región I")
    ax1.axvspan(r_minus, r_plus, alpha=0.16, label="Región II")
    ax1.set_xlim(1.8, 40)
    ax1.set_xlabel(r"$r/M$")
    ax1.set_ylabel(r"$V_{\rm eff},\,E^2$")
    ax1.set_title(r"Potencial efectivo, $L^2=16M^2$")
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(fontsize=7, ncol=2)

    ax2.plot(orbit.x/M, orbit.y/M, lw=2, label=rf"Trayectoria RK4, $r_0={r0:.1f}M$")
    draw_horizon(ax2)
    draw_circle(ax2, r_minus, label=rf"$r_-={r_minus:.1f}M$")
    draw_circle(ax2, r_plus, label=rf"$r_+={r_plus:.1f}M$")
    ax2.set_aspect("equal", "box")
    ax2.set_xlabel(r"$x/M$")
    ax2.set_ylabel(r"$y/M$")
    ax2.set_title("Trayectoria de caída directa")
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.32))
    savefig(fig, "figura3.png")


def figure4_bound_orbit():
    """Figura 4: órbita ligada con precesión, L^2=16M^2, r0=7M."""
    L2 = 16.0
    L = math.sqrt(L2)
    r0 = 7.0
    E = math.sqrt(v_eff(r0, L2, 1))
    orbit = integrate_orbit(E, L, 1, r0, 0.0, h=0.01, phi_goal=22*np.pi,
                            r_stop=2.02, max_steps=1_000_000)
    r_minus, r_plus = circular_radii_timelike(L2)
    r_plot = np.linspace(2.05, 40, 3000)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(r_plot/M, v_eff(r_plot, L2, 1), lw=2, label=r"$V_{\rm eff}$")
    ax1.axhline(E**2, ls="--", lw=1.3, label=rf"$E^2(r_0)={E**2:.3f}$")
    ax1.axvline(2, ls=":", lw=1.6, label=r"Horizonte $r_h=2M$")
    ax1.axvline(r_minus, ls="--", lw=1.5, label=rf"$r_-={r_minus:.1f}M$")
    ax1.axvline(r_plus, ls="--", lw=1.5, label=rf"$r_+={r_plus:.1f}M$")
    ax1.axvline(r0, ls=":", lw=1.5, label=rf"$r_0={r0:.1f}M$")
    ax1.axvspan(r_minus, r_plus, alpha=0.16, label="Región ligada")
    ax1.set_xlim(1.8, 40)
    ax1.set_xlabel(r"$r/M$")
    ax1.set_ylabel(r"$V_{\rm eff},\,E^2$")
    ax1.set_title(r"Potencial efectivo, $L^2=16M^2$")
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(fontsize=7, ncol=2)

    ax2.plot(orbit.x/M, orbit.y/M, lw=1.8, label=rf"Órbita ligada, $r_0={r0:.1f}M$")
    draw_horizon(ax2)
    draw_circle(ax2, r_minus, label=rf"$r_-={r_minus:.1f}M$")
    draw_circle(ax2, r_plus, label=rf"$r_+={r_plus:.1f}M$")
    ax2.set_aspect("equal", "box")
    ax2.set_xlabel(r"$x/M$")
    ax2.set_ylabel(r"$y/M$")
    ax2.set_title("Trayectoria ligada con precesión")
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.32))
    savefig(fig, "figura4.png")


def figure5_null_potential():
    """Figura 5: potencial efectivo nulo y esfera de fotones."""
    L2 = 36.0
    r_plot = np.linspace(2.05, 25, 2500)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(r_plot/M, v_eff(r_plot, L2, 0), lw=2, label=r"$V_{\rm eff}$")
    ax.axvline(2, ls="--", lw=1.5, label=r"$r_h=2M$ (horizonte)")
    ax.axvline(3, ls="--", lw=1.7, label=r"$r=3M$ (esfera de fotones)")
    ax.scatter([3], [v_eff(3, L2, 0)], s=55, zorder=5)
    ax.set_xlabel(r"$r/M$")
    ax.set_ylabel(r"$V_{\rm eff}$")
    ax.set_title(r"Potencial efectivo nulo, $\kappa=0$")
    ax.grid(True, ls=":", alpha=0.35)
    ax.legend(fontsize=8)
    savefig(fig, "figura5.png")


def integrate_null_from_infinity(b, y0, x0=20.0, h=0.01, max_steps=250_000):
    """Integra una trayectoria nula con E=1 y L=b."""
    E = 1.0
    L = b
    r0 = math.sqrt(x0**2 + y0**2)
    phi0 = math.atan2(y0, x0)
    radicand = E**2 - v_eff(r0, L**2, 0)
    if radicand <= 0:
        return np.array([]), np.array([]), "invalid_initial_data"
    ur0 = -math.sqrt(radicand)
    state = np.array([0.0, r0, ur0, phi0], dtype=float)
    xs, ys = [], []
    status = "max_steps"
    for step in range(max_steps):
        _, r, ur, phi = state
        if r <= 2.02:
            status = "captura"
            break
        if r > 35 and ur > 0 and step > 100:
            status = "dispersion"
            break
        if not np.all(np.isfinite(state)):
            status = "non_finite"
            break
        xs.append(r * math.cos(phi))
        ys.append(r * math.sin(phi))
        state = rk4_step(state, E, L, 0, h)
    return np.array(xs), np.array(ys), status


def figure6_null_classification():
    """Figura 6: clasificación de geodésicas nulas según b/bcrit."""
    bc = b_crit()
    cases = [
        (0.88*bc, "captura", r"$b<b_{\rm crit}$"),
        (1.00*bc, "crítica", r"$b=b_{\rm crit}$"),
        (1.25*bc, "dispersión", r"$b>b_{\rm crit}$"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, (b, title, label_b) in zip(axes, cases):
        x, y, _ = integrate_null_from_infinity(b=b, y0=b, x0=20, h=0.01)
        if len(x) > 0:
            ax.plot(x/M, y/M, lw=2, label=label_b)
        black_hole = plt.Circle((0, 0), 2, color="black")
        photon_sphere = plt.Circle((0, 0), 3, fill=False, ls="--", lw=1.3)
        ax.add_patch(black_hole)
        ax.add_patch(photon_sphere)
        ax.set_aspect("equal", "box")
        ax.set_xlabel(r"$x/M$")
        ax.set_ylabel(r"$y/M$")
        ax.set_title(title)
        ax.grid(True, ls=":", alpha=0.25)
        ax.legend(fontsize=8)
        if title == "dispersión":
            ax.set_xlim(-5, 22)
            ax.set_ylim(-20, 10)
        else:
            ax.set_xlim(-6.5, 6.5)
            ax.set_ylim(-6.5, 6.5)
    savefig(fig, "figura6.png")


# ============================================================
# Tablas de validación
# ============================================================

def generate_validation_tables():
    rows = []
    rows.append(("ISCO", r"$L^2=12M^2$", 6.0, 6.0, 0.0))
    r16 = circular_radii_timelike(16.0)
    rows.append((r"$r_-$", r"$L^2=16M^2$", r16[0], 4.0, abs(r16[0]-4.0)))
    rows.append((r"$r_+$", r"$L^2=16M^2$", r16[1], 12.0, abs(r16[1]-12.0)))
    r24 = circular_radii_timelike(24.0)
    rows.append((r"$r_-$", r"$L^2=24M^2$", r24[0], 3.5147186258, abs(r24[0]-3.5147186258)))
    rows.append((r"$r_+$", r"$L^2=24M^2$", r24[1], 20.4852813742, abs(r24[1]-20.4852813742)))
    rows.append((r"$r_{\rm ph}$", r"$\kappa=0$", photon_sphere_radius(), 3.0, 0.0))
    rows.append((r"$b_{\rm crit}$", r"$\kappa=0$", b_crit(), 5.1961524230, abs(b_crit()-5.1961524230)))

    conv_rows = []
    for h in [0.04, 0.02, 0.01, 0.005]:
        E = 1.0
        b = 7.0
        L = b * E
        r0 = 20.0
        ur0 = -math.sqrt(E**2 - v_eff(r0, L**2, 0))
        res = integrate_orbit(E, L, 0, r0, ur0, h=h, phi_goal=4*np.pi,
                              r_escape=30, max_steps=500_000)
        conv_rows.append((h, res.max_constraint_error, res.steps, res.status))

    txt = ["Validación de radios críticos en unidades M=1\n"]
    txt.append("Magnitud | Caso | Valor calculado | Valor analítico | Error absoluto\n")
    for mag, case, calc, analytic, err in rows:
        txt.append(f"{mag} | {case} | {calc:.10f} | {analytic:.10f} | {err:.3e}\n")
    txt.append("\nConservación de la relación radial para una geodésica nula de dispersión\n")
    txt.append("h | Delta_max | pasos | estado\n")
    for h, err, steps, status in conv_rows:
        txt.append(f"{h:.3f} | {err:.3e} | {steps} | {status}\n")
    (TAB_DIR / "tabla_validacion_RMF.txt").write_text("".join(txt), encoding="utf-8")

    latex = r"""% Tabla generada automáticamente por schwarzschild_geodesics_RMF_clean.py

\begin{table}[H]
\centering
\caption{Validaci\'on de radios cr\'iticos y par\'ametros caracter\'isticos en unidades $M=1$.}
\label{tab:validacion_radios}
\scriptsize
\renewcommand{\arraystretch}{1.08}
\begin{tabular}{c c c c}
\hline
Magnitud & Caso & Valor & Error abs. \\
\hline
"""
    for mag, case, calc, analytic, err in rows:
        latex += f"{mag} & {case} & {calc:.6f} & ${err:.1e}$ \\\\n"
    latex += r"""\hline
\end{tabular}
\normalsize
\end{table}

\begin{table}[H]
\centering
\caption{Conservaci\'on de la relaci\'on radial $\Delta=\max|u_r^2+V_{\mathrm{eff}}(r)-E^2|$ para una geod\'esica nula de dispersi\'on con $E=1$, $b=7M$ y $r_0=20M$.}
\label{tab:convergencia_rk4}
\scriptsize
\renewcommand{\arraystretch}{1.08}
\begin{tabular}{c c c c}
\hline
$h$ & $\Delta_{\max}$ & Pasos & Estado \\
\hline
"""
    for h, err, steps, status in conv_rows:
        latex += f"{h:.3f} & ${err:.3e}$ & {steps} & {status} \\\\n"
    latex += r"""\hline
\end{tabular}
\normalsize
\end{table}
"""
    (TAB_DIR / "tabla_validacion_RMF.tex").write_text(latex, encoding="utf-8")


def generate_readme_and_requirements():
    readme = """# Schwarzschild Geodesics RMF

Código reproducible para generar las figuras y tablas del manuscrito preparado para Revista Mexicana de Física.

## Ejecución

```bash
pip install -r requirements.txt
python schwarzschild_geodesics_RMF_clean.py
```

Los resultados se guardan en:

```text
output_rmf/figures/
output_rmf/tables/
```

## Autores

- Ricardo Angelo Ballon Tito
- Rolando Moisés Perca Gonzales
"""
    (OUTPUT_DIR / "README.md").write_text(readme, encoding="utf-8")
    (OUTPUT_DIR / "requirements.txt").write_text("numpy\nmatplotlib\n", encoding="utf-8")


def main():
    print("Generando figuras y tablas para RMF...")
    figure1_isco()
    figure2_l16_potential()
    figure3_direct_infall()
    figure4_bound_orbit()
    figure5_null_potential()
    figure6_null_classification()
    generate_validation_tables()
    generate_readme_and_requirements()
    print("Listo.")
    print(f"Figuras: {FIG_DIR.resolve()}")
    print(f"Tablas:  {TAB_DIR.resolve()}")


if __name__ == "__main__":
    main()
