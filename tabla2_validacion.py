import os
import numpy as np

M = 1.0

def f(r):
    return 1 - 2*M/r

def dV_dr(r, L2, kappa):
    return (2*M/r**2)*(kappa + L2/r**2) + f(r)*(-2*L2/r**3)

def d2V_dr2(r, L2, kappa):
    return -4*M*kappa/r**3 + 6*L2/r**4 - 24*M*L2/r**5

def biseccion(fun, a, b, tol=1e-12, max_iter=500):
    fa = fun(a)
    fb = fun(b)

    if fa*fb > 0:
        raise ValueError('El intervalo no contiene una raíz.')

    for _ in range(max_iter):
        c = 0.5*(a + b)
        fc = fun(c)

        if abs(b - a)/2 < tol:
            return c

        if fa*fc <= 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

    return 0.5*(a + b)

def radios_analiticos(L2):
    raiz = np.sqrt(1 - 12*M**2/L2)
    r_minus = (L2/(2*M))*(1 - raiz)
    r_plus = (L2/(2*M))*(1 + raiz)
    return r_minus, r_plus

def L2_circular(r):
    return M*r**2/(r - 3*M)

r_isco_num = biseccion(lambda r: d2V_dr2(r, L2_circular(r), 1), 5.0, 7.0)

r16_minus_num = biseccion(lambda r: dV_dr(r, 16.0, 1), 2.01, 5.99)
r16_plus_num = biseccion(lambda r: dV_dr(r, 16.0, 1), 6.01, 30.0)
r16_minus_an, r16_plus_an = radios_analiticos(16.0)

r24_minus_num = biseccion(lambda r: dV_dr(r, 24.0, 1), 2.01, 5.99)
r24_plus_num = biseccion(lambda r: dV_dr(r, 24.0, 1), 6.01, 40.0)
r24_minus_an, r24_plus_an = radios_analiticos(24.0)

r_ph_num = biseccion(lambda r: dV_dr(r, 36.0, 0), 2.01, 10.0)
b_crit_num = r_ph_num/np.sqrt(f(r_ph_num))

datos = [
    ('ISCO', 6.0, r_isco_num),
    ('r_-  L2=16', r16_minus_an, r16_minus_num),
    ('r_+  L2=16', r16_plus_an, r16_plus_num),
    ('r_-  L2=24', r24_minus_an, r24_minus_num),
    ('r_+  L2=24', r24_plus_an, r24_plus_num),
    ('r_ph', 3.0, r_ph_num),
    ('b_crit', 3*np.sqrt(3), b_crit_num),
]

os.makedirs('output/tables', exist_ok=True)

lineas = ['Magnitud | Analítico | Numérico | Error absoluto\n']
for nombre, analitico, numerico in datos:
    error = abs(numerico - analitico)
    lineas.append(f'{nombre} | {analitico:.10f} | {numerico:.10f} | {error:.3e}\n')

with open('output/tables/tabla2_validacion.txt', 'w', encoding='utf-8') as archivo:
    archivo.writelines(lineas)

for linea in lineas:
    print(linea.strip())
