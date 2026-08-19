import os
import numpy as np
import matplotlib.pyplot as plt

# Parámetros
M = 1.0
L = 6.0
r_h = 2.0*M
r_ph = 3.0*M

def V_std(r, M, L):
    return (1.0 - 2.0*M/r)*(L**2)/(r**2)

r = np.linspace(2.01*M, 25.0*M, 4000)
V = V_std(r, M, L)

plt.figure(figsize=(8, 4.8))
plt.plot(r/M, V, color='blue', lw=2)
plt.axvline(r_ph/M, linestyle='--', color='orange',
            label=r'$r=3M$ (esfera de fotones)')
plt.axvline(r_h/M, linestyle='--', color='black',
            label=r'$r_h=2M$ (horizonte)')

current_ticks = plt.gca().get_xticks()
ticks = np.unique(np.concatenate([current_ticks, [r_h/M, r_ph/M]]))
labels = []
for t in ticks:
    if np.isclose(t, r_h/M):
        labels.append('2')
    elif np.isclose(t, r_ph/M):
        labels.append('3')
    else:
        labels.append(f'{t:g}')

plt.xticks(ticks, labels)
plt.xlabel(r'$r/M$')
plt.ylabel(r'$V_{\rm eff}$')
plt.title(r'Potencial efectivo nulo, $\kappa=0$')
plt.legend()
plt.grid(True, ls=':', alpha=0.25)
plt.tight_layout()

os.makedirs('output/figures', exist_ok=True)
plt.savefig('output/figures/figura5.png', dpi=300, bbox_inches='tight')
plt.close()
