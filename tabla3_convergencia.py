import os
import numpy as np

M = 1.0
E = 1.0
b = 7.0
L = b
r0 = 20.0

def f(r):
    return 1 - 2*M/r

def V_eff(r):
    return f(r)*(L**2/r**2)

def dV_dr(r):
    return (2*M/r**2)*(L**2/r**2) + f(r)*(-2*L**2/r**3)

def rhs(s):
    t, r, pr, phi = s
    return np.array([E/f(r), pr, -0.5*dV_dr(r), L/r**2])

def rk4_step(s, h):
    k1 = rhs(s)
    k2 = rhs(s + 0.5*h*k1)
    k3 = rhs(s + 0.5*h*k2)
    k4 = rhs(s + h*k3)
    return s + h*(k1 + 2*k2 + 2*k3 + k4)/6

def error_restriccion(r, pr):
    return abs(pr**2 + V_eff(r) - E**2)

def integrar(h):
    pr0 = -np.sqrt(E**2 - V_eff(r0))
    st = np.array([0.0, r0, pr0, 0.0])

    error_max = 0.0
    pasos = 0
    estado = 'máximo de pasos'

    for i in range(500000):
        t, r, pr, phi = st

        error = error_restriccion(r, pr)
        if error > error_max:
            error_max = error

        if r <= 2.02*M:
            estado = 'captura'
            break

        if i > 20 and r >= 30*M and pr > 0:
            estado = 'dispersión'
            break

        if phi >= 4*np.pi:
            estado = 'ángulo máximo'
            break

        st = rk4_step(st, h)
        pasos += 1

    return error_max, pasos, estado

os.makedirs('output/tables', exist_ok=True)
hs = [0.04, 0.02, 0.01, 0.005]

lineas = ['h | Delta_max | pasos | estado\n']
for h in hs:
    error, pasos, estado = integrar(h)
    lineas.append(f'{h:.3f} | {error:.3e} | {pasos} | {estado}\n')

with open('output/tables/tabla3_convergencia.txt', 'w', encoding='utf-8') as archivo:
    archivo.writelines(lineas)

for linea in lineas:
    print(linea.strip())
