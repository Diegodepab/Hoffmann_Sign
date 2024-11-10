#!/bin/bash

# Definir las carpetas donde estarán las dependencias de python y de R
export R_LIB: ./R_packages
mkdir_p$R_LIB
R.script_e "install.packages(Pack)"

export PYTHON_LIB=./py_packages
pip install -r requirements.txt

# Limpieza de cache
python -m pip cache purge

# Instalar dependencias de R del archivo Rreqs.txt en la carpeta R
Rscript -e 'install.packages(readLines("Rreqs.txt"), lib="./R_packages")'

# Imprimir mensaje de finalización
echo "Instalación de dependencias completada con éxito.

