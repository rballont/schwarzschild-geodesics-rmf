import subprocess
import sys

scripts = [
    'figura1_isco.py',
    'figura2_potencial_L16.py',
    'figura3_caida_r35.py',
    'figura4_orbita_ligada.py',
    'figura5_potencial_nulo.py',
    'figura6_captura_critica_dispersion.py',
    'tabla2_validacion.py',
    'tabla3_convergencia.py',
]

for script in scripts:
    print('\nEjecutando:', script)
    subprocess.run([sys.executable, script], check=True)

print('\nTodos los resultados fueron generados.')
