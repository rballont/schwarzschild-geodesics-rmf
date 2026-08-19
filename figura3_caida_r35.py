import os
import numpy as np
import matplotlib.pyplot as plt

# ========= 1. parámetros físicos =========
M = 1.0
L2 = 16.0
L = np.sqrt(L2)

def f(r):
    return 1 - 2*M/r

def V_eff(r):
    return f(r)*(1 + L2/r**2)

def dV_dr(r):
    return 2*M/r**2 + (6*M*L2)/r**4 - (2*L2)/r**3

# ========= 2. radios circulares =========
disc_term = 1 - 12*M**2/L2
sqrt_disc = np.sqrt(disc_term)
r_minus = (L2/(2*M))*(1 - sqrt_disc)
r_plus = (L2/(2*M))*(1 + sqrt_disc)

# ========= 3. Runge-Kutta =========
def rhs(s, E):
    t, r, pr, phi = s
    return np.array([E/f(r), pr, -0.5*dV_dr(r), L/r**2])

def rk4_step(s, E, h=0.002):
    k1 = rhs(s, E)
    k2 = rhs(s + 0.5*h*k1, E)
    k3 = rhs(s + 0.5*h*k2, E)
    k4 = rhs(s + h*k3, E)
    return s + h*(k1 + 2*k2 + 2*k3 + k4)/6

r0 = 3.5*M
E = np.sqrt(V_eff(r0))
st = np.array([0.0, r0, 0.0, 0.0])
xs, ys = [], []

for _ in range(500000):
    t, r, pr, phi = st
    xs.append(r*np.cos(phi))
    ys.append(r*np.sin(phi))
    if r <= 2.02*M:
        break
    if phi >= 6*np.pi:
        break
    st = rk4_step(st, E, h=0.002)

# ========= 4. figura =========
r_mesh = np.linspace(2.01*M, 40*M, 2000)
theta = np.linspace(0, 2*np.pi, 400)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

ax1.plot(r_mesh/M, V_eff(r_mesh), lw=2, color='tab:orange',
         label=r'$V_{\mathrm{eff}}$')
ax1.axvspan(2, r_minus/M, color='#f4b6c2', alpha=0.30,
            label=fr'Región I: $2M<r<{r_minus/M:.1f}M$')
ax1.axvspan(r_minus/M, r_plus/M, color='#00a86b', alpha=0.20,
            label=fr'Región II: ${r_minus/M:.1f}M\leq r\leq {r_plus/M:.1f}M$')
ax1.axvline(r_minus/M, ls='--', color='#006400', label=fr'$r_-={r_minus/M:.1f}M$')
ax1.axvline(r_plus/M, ls='--', color='#800080', label=fr'$r_+={r_plus/M:.1f}M$')
ax1.axhline(E**2, ls='--', color='tab:red', label=fr'$E^2(r_0)={E**2:.3f}$')
ax1.scatter(r0/M, E**2, color='tab:red', s=60, zorder=5,
            label=fr'$r_0={r0/M:.1f}M$')
ax1.axvline(2, ls=':', lw=1.8, color='black', label=r'Horizonte $r_h=2M$')
ax1.set_xlim(1.8, 40)
ax1.set_xlabel(r'$r/M$')
ax1.set_ylabel(r'$V_{\mathrm{eff}},\,E^2$')
ax1.set_title(rf'Potencial efectivo ($L^2={L2:.0f}$)')
ax1.grid(True, ls=':', alpha=0.25)
ax1.legend(fontsize=8, ncol=2)

ax2.plot(np.array(xs)/M, np.array(ys)/M, color='tab:red', lw=1.8,
         label=fr'Trayectoria $r_0={r0/M:.1f}M$')
ax2.plot(2*np.cos(theta), 2*np.sin(theta), 'k--', label=r'Horizonte $r=2M$')
ax2.set_aspect('equal', 'box')
ax2.set_xlabel(r'$x/M$')
ax2.set_ylabel(r'$y/M$')
ax2.set_title('Trayectoria de caída')
ax2.grid(True, ls=':', alpha=0.25)
ax2.legend(fontsize=8)

plt.tight_layout()
os.makedirs('output/figures', exist_ok=True)
plt.savefig('output/figures/figura3.png', dpi=300, bbox_inches='tight')
plt.close()
