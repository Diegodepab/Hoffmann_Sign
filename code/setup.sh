#!/bin/bash

# En caso de que ocurra un error, detener la ejecución del script
set -e

# Definir las carpetas donde estarán las dependencias de Python y de R
export R_LIB="./R_packages"
mkdir -p $R_LIB

export PYTHON_LIB="./py_packages"
mkdir -p $PYTHON_LIB

# Instalar dependencias de Python en el directorio personalizado
pip install --prefix=$PYTHON_LIB -r requirements.txt

# Limpieza de caché
python -m pip cache purge

# Instalar dependencias de R del archivo Rreqs.txt en la carpeta R
Rscript -e 'install.packages(readLines("Rreqs.txt"), lib="./R_packages")'

# Imprimir mensaje de finalización
echo "Instalación de dependencias completada con éxito."
