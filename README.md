# Schwarzschild Geodesics — RMF

Códigos reproducibles del manuscrito:

**“Clasificación analítico-computacional de geodésicas tipo tiempo y nulas en el espacio-tiempo de Schwarzschild”**

Autores:
- Ricardo Angelo Ballon Tito
- Rolando Moisés Perca Gonzales

## Relación con la tesis

Los códigos de las figuras provienen de los programas desarrollados por Ricardo Angelo Ballon Tito en el **Anexo D** de su tesis de licenciatura. La estructura computacional corresponde a la empleada en ese trabajo: scripts autocontenidos, integración RK4 escrita explícitamente y dependencias centradas en `numpy` y `matplotlib`.

| Manuscrito RMF | Código | Base en la tesis |
|---|---|---|
| Figura 1 | `figura1_isco.py` | Figura 8 — ISCO, L²=12M² |
| Figura 2 | `figura2_potencial_L16.py` | Figura 9 — potencial, L²=16M² |
| Figura 3 | `figura3_caida_r35.py` | Figura 14 — región I |
| Figura 4 | `figura4_orbita_ligada.py` | Figura 15 — región II |
| Figura 5 | `figura5_potencial_nulo.py` | Figura 17 — potencial nulo |
| Figura 6 | `figura6_captura_critica_dispersion.py` | Figura 19 — geodésicas nulas |

La **Figura 6** deriva del código de la Figura 19 del Anexo D. En el manuscrito RMF se emplean las condiciones E=1, r0=30M y b/bcrit=0.9, 1.0, 1.2. La convención gráfica corresponde a la usada en la tesis: región del agujero negro en negro y trayectorias de los fotones en amarillo.

Las Tablas II y III corresponden a comprobaciones numéricas del manuscrito: concordancia analítica-numérica y conservación de la restricción radial.

## Dependencias

```bash
pip install -r requirements.txt
```

Se requieren:

```text
numpy
matplotlib
```

## Ejecución por figura

```bash
python figura1_isco.py
python figura2_potencial_L16.py
python figura3_caida_r35.py
python figura4_orbita_ligada.py
python figura5_potencial_nulo.py
python figura6_captura_critica_dispersion.py
```

Para las tablas:

```bash
python tabla2_validacion.py
python tabla3_convergencia.py
```

Para reproducir todo:

```bash
python run_all.py
```

Las figuras se guardan en `output/figures/` y las tablas en `output/tables/`.

## Unidades

Se emplean unidades geométricas: G=c=1 y M=1.

## Licencia

MIT.
