#!/bin/bash

# En caso de que ocurra un error, detener la ejecución del script
set -e

# Definir las carpetas donde estarán las dependencias de Python y de R
export R_LIB="./R_packages"
mkdir -p $R_LIB

python_libs="./py_packages"
export PYTHONPATH=$python_libs:$PYTHONPATH
mkdir -p $python_libs

# Verificar si R está instalado
if ! command -v Rscript &>/dev/null; then
  echo "R no está instalado. Por favor, instala R antes de continuar."
  exit 1
fi

# Verificar si Python está instalado
if ! command -v python &>/dev/null; then
  echo "Python no está instalado. Por favor, instala Python antes de continuar."
  exit 1
fi

# Instalar dependencias de Python en el directorio personalizado
pip install -r requirements.txt

# Limpieza de caché de Python
#python -m pip cache purge

# Instalar dependencias de R desde el archivo Rreqs.txt en la carpeta R
export R_LIBS_USER="./R_packages"
Rscript -e 'install.packages(readLines("Rreqs.txt"), lib='$R_LIBS_USER', repos="https://cloud.r-project.org")'
# Instalar dependencias de Bioconductor desde el archivo BioconductorReqs.txt en la carpeta R
Rscript -e 'BiocManager::install(readLines("BioconductorReqs.txt"), lib='$R_LIBS_USER')'
# Limpiar caché de R (opcional, puede ayudar a prevenir errores de instalación)
#Rscript -e 'unlink(file.path(Sys.getenv("R_LIBS_USER"), "00LOCK"), recursive = TRUE)'

# Imprimir mensaje de finalización
echo "Instalación de dependencias completada con éxito."


