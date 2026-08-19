import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------ Parámetros globales ------------------
M = 1.0
E = 1.0
r_h = 2.0*M
r_ph = 3.0*M
b_crit = 3.0*np.sqrt(3.0)*M

# En el artículo los tres casos parten desde el mismo radio.
r0 = 30.0*M
h = 0.002

# ------------------ Utilidades ------------------
def f_metric(r, M):
    return 1.0 - 2.0*M/r

def derivs(lam, y, M, E, L):
    # y = [t, r, phi, pr], con pr = r'
    t, r, phi, pr = y
    if r <= 2.02*M:
        r = 2.02*M

    f = f_metric(r, M)
    t_dot = E/f
    phi_dot = L/(r**2)

    # Forma radial equivalente obtenida del potencial efectivo nulo:
    # r'' = -1/2 dV_eff/dr = L^2/r^3 - 3 M L^2/r^4
    r_ddot = L**2/r**3 - 3*M*L**2/r**4

    return np.array([t_dot, pr, phi_dot, r_ddot])

def rk4(fun, lam, y, h, *args):
    k1 = fun(lam, y, *args)
    k2 = fun(lam + 0.5*h, y + 0.5*h*k1, *args)
    k3 = fun(lam + 0.5*h, y + 0.5*h*k2, *args)
    k4 = fun(lam + h, y + h*k3, *args)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def integrate_case(factor_b, N=2000000):
    b = factor_b*b_crit
    L = b*E

    f0 = f_metric(r0, M)
    pr2 = E**2 - f0*(L**2)/(r0**2)
    if pr2 <= 0:
        raise ValueError('Condiciones iniciales inconsistentes.')

    pr0 = -np.sqrt(pr2)
    y = np.array([0.0, r0, 0.0, pr0])
    lam = 0.0
    Y = [y.copy()]

    for i in range(N):
        r = y[1]
        pr = y[3]

        if factor_b < 1.0 and r <= r_h*(1 + 1e-3):
            break

        if np.isclose(factor_b, 1.0) and r <= 3.00005*M and abs(pr) < 5e-4:
            break

        if factor_b > 1.0 and i > 500 and r >= r0 and pr > 0:
            break

        if y[2] > 20*np.pi:
            break

        y = rk4(derivs, lam, y, h, M, E, L)
        Y.append(y.copy())
        lam += h

    return np.array(Y)

def to_cartesian(rs, phis):
    xs = rs*np.cos(phis)
    ys = rs*np.sin(phis)
    return xs, ys

# ------------------ Casos del manuscrito ------------------
factores = [0.9, 1.0, 1.2]
titulos = [
    r'(a) Captura: $b=0.9\,b_{\rm crit}$',
    r'(b) Trayectoria crítica: $b=b_{\rm crit}$',
    r'(c) Dispersión: $b=1.2\,b_{\rm crit}$'
]

trayectorias = []
for factor in factores:
    Y = integrate_case(factor)
    t, r, phi, pr = Y.T
    x, y = to_cartesian(r, phi)
    trayectorias.append((x, y))

# ------------------ Gráficos ------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
yellow = '#FFD400'
theta = np.linspace(0, 2*np.pi, 500)

for ax, (x, y), titulo in zip(axes, trayectorias, titulos):
    bh = plt.Circle((0, 0), r_h/M, facecolor='black',
                    edgecolor='black', linewidth=1.4, zorder=1)
    ax.add_patch(bh)

    ax.plot((r_ph/M)*np.cos(theta), (r_ph/M)*np.sin(theta),
            '--', color=yellow, lw=1.4, label=r'Esfera de fotones $r=3M$')

    ax.plot(x/M, y/M, color=yellow, lw=2.2, zorder=2)

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(r'$x/M$')
    ax.set_ylabel(r'$y/M$')
    ax.set_title(titulo)
    ax.grid(True, ls=':', alpha=0.25)

axes[0].set_xlim(-6, 31)
axes[0].set_ylim(-8, 12)
axes[1].set_xlim(-10, 31)
axes[1].set_ylim(-12, 12)
axes[2].set_xlim(-8, 31)
axes[2].set_ylim(-18, 14)

plt.tight_layout()
os.makedirs('output/figures', exist_ok=True)
plt.savefig('output/figures/figura6.png', dpi=300, bbox_inches='tight')
plt.close()

print('b_crit =', b_crit)
for factor, (x, y) in zip(factores, trayectorias):
    print(f'b/b_crit = {factor:.1f}   puntos = {len(x)}')
