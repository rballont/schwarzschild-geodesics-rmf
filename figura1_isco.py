import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------- parámetros físicos ----------------------
M = 1.0
L2 = 12.0
L = np.sqrt(L2)

def f(r):
    return 1 - 2*M/r

def V_eff(r):
    return f(r)*(1 + L2/r**2)

def dV_dr(r):
    return (2*M/r**2)*(1 + L2/r**2) + f(r)*(-2*L2/r**3)

# ------------------- condiciones iniciales -------------------
r0 = 6.0*M
E = np.sqrt(V_eff(r0))
state = np.array([0.0, r0, 0.0, 0.0])  # (t, r, p_r, phi)

# ------------------- integrador RK-4 -------------------------
def rhs(s):
    t, r, pr, phi = s
    return np.array([E/f(r), pr, -0.5*dV_dr(r), L/r**2])

def rk4(s, h):
    k1 = rhs(s)
    k2 = rhs(s + 0.5*h*k1)
    k3 = rhs(s + 0.5*h*k2)
    k4 = rhs(s + h*k3)
    return s + h*(k1 + 2*k2 + 2*k3 + k4)/6

h = 0.01  # valor usado en el manuscrito RMF
phi_max = 20*np.pi
xs, ys = [], []

while state[3] < phi_max:
    t, r, pr, phi = state
    xs.append(r*np.cos(phi))
    ys.append(r*np.sin(phi))
    state = rk4(state, h)

# ------------------- figura con dos paneles ------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

r_plot = np.linspace(2.2*M, 20*M, 1400)
ax1.plot(r_plot/M, V_eff(r_plot), lw=2, label=r'$V_{\mathrm{eff}}$')
ax1.axhline(E**2, ls='--', color='tab:red', label=fr'$E^2={E**2:.4f}$')
ax1.scatter(r0/M, E**2, color='tab:red', s=70, zorder=5, label=r'$r_0=6M$')
ax1.axvline(6, ls='--', lw=2, color='gold', label=r'ISCO $r=6M$')
ax1.axvline(2, ls='--', lw=1.8, color='0.25', label=r'Horizonte $r=2M$')
ax1.set_xlim(1.8, 20)
ax1.set_xticks(np.arange(2, 21, 2))
ax1.set_xlabel(r'$r/M$')
ax1.set_ylabel(r'$V_{\mathrm{eff}},\,E^2$')
ax1.set_title('Potencial efectivo')
ax1.legend(fontsize='small')
ax1.grid(True, ls=':', alpha=0.25)

th = np.linspace(0, 2*np.pi, 400)
ax2.plot(xs, ys, color='tab:red', lw=2, label='RK-4')
ax2.plot(2*M*np.cos(th), 2*M*np.sin(th), 'm--', label=r'Horizonte $r=2M$')
ax2.set_aspect('equal', 'box')
ax2.set_xlabel(r'$x/M$')
ax2.set_ylabel(r'$y/M$')
ax2.set_title('Órbita circular ISCO')
ax2.legend(fontsize='small')
ax2.grid(True, ls=':', alpha=0.25)

plt.tight_layout()
os.makedirs('output/figures', exist_ok=True)
plt.savefig('output/figures/figura1.png', dpi=300, bbox_inches='tight')
plt.close()
