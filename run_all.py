#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Ejecuta todos los códigos asociados a las figuras y tablas del manuscrito."""

import subprocess
import sys

SCRIPTS = [
    "figura1_isco.py",
    "figura2_potencial_L16.py",
    "figura3_caida_r35.py",
    "figura4_orbita_ligada.py",
    "figura5_potencial_nulo.py",
    "figura6_captura_critica_dispersion.py",
    "tabla2_validacion.py",
    "tabla3_convergencia.py",
]


def main():
    for script in SCRIPTS:
        print(f"Ejecutando {script}...")
        subprocess.run([sys.executable, script], check=True)
    print("Reproducción completa finalizada.")


if __name__ == "__main__":
    main()
