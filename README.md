# Schwarzschild Geodesics — RMF

Código reproducible asociado al manuscrito **“Clasificación analítico-computacional de geodésicas tipo tiempo y nulas en el espacio-tiempo de Schwarzschild”**, preparado para la Revista Mexicana de Física.

## Autores

- Ricardo Angelo Ballon Tito — Universidad Nacional de San Agustín de Arequipa — ORCID: 0000-0002-5438-7699
- Rolando Moisés Perca Gonzales — Departamento Académico de Física, Universidad Nacional de San Agustín de Arequipa — ORCID: 0000-0001-5734-7310

## Requisitos

```bash
pip install -r requirements.txt
```

El código usa unidades geométricas `G = c = 1` y en las simulaciones del manuscrito se toma `M = 1`.

## Reproducción por figura

Cada figura del artículo tiene un script independiente:

| Resultado del manuscrito | Código | Salida |
|---|---|---|
| Figura 1 — ISCO | `figura1_isco.py` | `output/figures/figura1.png` |
| Figura 2 — Potencial, `L²=16M²` | `figura2_potencial_L16.py` | `output/figures/figura2.png` |
| Figura 3 — Caída desde `r0=3.5M` | `figura3_caida_r35.py` | `output/figures/figura3.png` |
| Figura 4 — Órbita ligada | `figura4_orbita_ligada.py` | `output/figures/figura4.png` |
| Figura 5 — Potencial nulo | `figura5_potencial_nulo.py` | `output/figures/figura5.png` |
| Figura 6 — Captura, crítica y dispersión | `figura6_captura_critica_dispersion.py` | `output/figures/figura6.png` |
| Tabla II — Validación de radios | `tabla2_validacion.py` | `output/tables/tabla2_validacion.txt` |
| Tabla III — Convergencia RK4 | `tabla3_convergencia.py` | `output/tables/tabla3_convergencia.txt` |

Las funciones que comparten los distintos casos —métrica de Schwarzschild, potencial efectivo, derivadas, integración RK4 y búsqueda de raíces— se encuentran en `schwarzschild_utils.py`.

## Ejecutar una sola figura

Ejemplo:

```bash
python figura1_isco.py
```

Para la Figura 6:

```bash
python figura6_captura_critica_dispersion.py
```

En esta última se usan las condiciones del manuscrito:

- `E = 1`
- `r0 = 30M`
- `h = 0.002`
- `b/bcrit = 0.9, 1.0, 1.2`

La trayectoria de los fotones se representa en amarillo/dorado y el agujero negro en negro, manteniendo el estilo gráfico usado en el manuscrito.

## Ejecutar todos los resultados

```bash
python run_all.py
```

Las salidas se crean en:

```text
output/
├── figures/
└── tables/
```

## Validación numérica

`tabla2_validacion.py` localiza numéricamente los ceros de las condiciones diferenciales del potencial mediante bisección con tolerancia `1e-12`. Los valores se comparan con las expresiones analíticas usadas en el manuscrito.

`tabla3_convergencia.py` evalúa

```text
Delta = max |u_r^2 + V_eff(r) - E^2|
```

para una geodésica nula de dispersión con `E=1`, `b=7M` y `r0=20M`, usando distintos tamaños de paso RK4.

## Archivo histórico

`schwarzschild_geodesics_RMF_clean.py` corresponde a una versión previa en la que todas las figuras y tablas se generaban desde un único archivo. Se conserva como referencia histórica; para reproducir la versión actual del manuscrito deben usarse los scripts individuales indicados arriba.

## Licencia

MIT License.
