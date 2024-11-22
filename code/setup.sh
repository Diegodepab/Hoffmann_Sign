#!/bin/bash

# En caso de que ocurra un error, detener la ejecución del script
set -e

# Definir las carpetas donde estarán las dependencias de Python y de R
export R_LIB="./R_packages"
mkdir -p $R_LIB

export PYTHON_LIB="./py_packages"
mkdir -p $PYTHON_LIB

# Verificar si R está instalado
if ! command -v Rscript &>/dev/null; then
  echo "R no está instalado. Por favor, instala R antes de continuar."
  exit 1
fi

# Verificar si Python está instalado
if ! command -v python3 &>/dev/null; then
  echo "Python no está instalado. Por favor, instala Python antes de continuar."
  exit 1
fi

# Instalar dependencias de Python en el directorio personalizado
pip install --prefix=$PYTHON_LIB -r requirements.txt

# Limpieza de caché de Python
python -m pip cache purge

# Instalar dependencias de R desde el archivo Rreqs.txt en la carpeta R
Rscript -e 'install.packages(readLines("Rreqs.txt"), lib="./R_packages")'

# Instalar librerías necesarias de R si no están incluidas en Rreqs.txt
Rscript -e 'if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph", lib="./R_packages")'
Rscript -e 'if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("clusterProfiler", lib="./R_packages")'
Rscript -e 'if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("org.Hs.eg.db", lib="./R_packages")'

# Limpiar caché de R (opcional, puede ayudar a prevenir errores de instalación)
Rscript -e 'unlink(file.path(Sys.getenv("R_LIBS_USER"), "00LOCK"), recursive = TRUE)'

# Imprimir mensaje de finalización
echo "Instalación de dependencias completada con éxito."

