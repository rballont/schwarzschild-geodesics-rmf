import os
import numpy as np
import matplotlib.pyplot as plt

# 1. Parámetros -------------------------------------------------
M = 1.0
L_sq = 16.0

# 2. Potencial y derivada --------------------------------------
def V_eff(r):
    return (1 - 2*M/r)*(1 + L_sq/r**2)

def dV_dr(r):
    return 2*M/r**2 - 2*L_sq/r**3 + 6*M*L_sq/r**4

# 3. Extremos ---------------------------------------------------
disc = L_sq**2 - 12*M**2*L_sq
r_minus = (L_sq - np.sqrt(disc))/(2*M)
r_plus = (L_sq + np.sqrt(disc))/(2*M)
V_max = V_eff(r_minus)
V_min = V_eff(r_plus)

# 4. Gráfica ----------------------------------------------------
r = np.linspace(2.05*M, 30*M, 3000)
fig, ax = plt.subplots(figsize=(9, 4.8))

ax.plot(r/M, V_eff(r), lw=2, label=r'$V_{\mathrm{eff}}(r)$')
ax.axhline(1, ls=':', lw=1, color='black', label=r'$E^2=1$')

ax.scatter([r_minus/M, r_plus/M], [V_max, V_min],
           color='tab:orange', zorder=5)
ax.annotate(rf'$r_-={r_minus/M:.1f}M$', (r_minus/M, V_max),
            xytext=(6, 8), textcoords='offset points', fontsize=9)
ax.annotate(rf'$r_+={r_plus/M:.1f}M$', (r_plus/M, V_min),
            xytext=(6, -15), textcoords='offset points', fontsize=9)

rosa_palo = '#f4b6c2'
verde_jade = '#00a86b'
ax.axvspan(2, r_minus/M, color=rosa_palo, alpha=0.15,
           label='Región I: caída directa')
ax.axvspan(r_minus/M, r_plus/M, color=verde_jade, alpha=0.12,
           label='Región II: órbitas posibles')

ax.axvline(2, ls=':', lw=1.8, color='black', label=r'Horizonte $r_h=2M$')
ax.axvline(r_minus/M, ls='--', color=rosa_palo, label=r'$r_-$ (máx.)')
ax.axvline(r_plus/M, ls='--', color=verde_jade, label=r'$r_+$ (mín.)')

ax.set_xlabel(r'$r/M$')
ax.set_ylabel(r'$V_{\mathrm{eff}}(r)$')
ax.set_title(rf'Potencial efectivo tipo tiempo ($M=1$, $L^2={L_sq:.0f}$)')
ax.set_xlim(1.8, 30)
ax.set_xticks(np.arange(2, 31, 2))
ax.set_ylim(0, 1.15)
ax.grid(alpha=0.3, linestyle=':')
ax.legend(frameon=False, fontsize=8, loc='upper center',
          bbox_to_anchor=(0.5, -0.17), ncol=3)

fig.tight_layout()
os.makedirs('output/figures', exist_ok=True)
plt.savefig('output/figures/figura2.png', dpi=300, bbox_inches='tight')
plt.close()

print(f'r_- (máximo inestable): {r_minus:.6f} M   V_max = {V_max:.6f}')
print(f'r_+ (mínimo estable) : {r_plus:.6f} M   V_min = {V_min:.6f}')
