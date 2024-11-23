#!/bin/bash

# Definir las carpetas donde estarán las dependencias de Python y de R
export R_LIB="./R_packages"
export PYTHONPATH="./py_packages/local/lib/python3.11/dist-packages:$PYTHONPATH"

# Crear directorios necesarios
mkdir -p data results images logs

# Descargar los datos de genes de HPO
echo "Descargando datos de genes de HPO..."
wget -q https://ontology.jax.org/api/network/annotation/HP:0031993/download/gene -O data/genes.tsv

# Convertir genes a string ID (descomentar cuando el script esté listo)
echo "Convirtiendo genes a STRING ID..."
python genes2string.py --genes_iniciales data/genes.tsv --genes_mapeados data/genes_string.tsv

# Descargar la red de String
echo "Descargando red STRING..."
wget -q https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz -O data/network.txt.gz

# Descomprimir el archivo descargado
echo "Descomprimiendo red STRING..."
gunzip -f data/network.txt.gz

# Propagación de genes con DIAMOnD
echo "Propagando genes usando DIAMOnD..."
python DIAMOnD.py data/genes_string.tsv data/network.txt data/red_propagada.txt

# Descargar interacciones usando el script string_interactions.py
echo "Obteniendo interacciones de STRING..."
python string_interactions.py --input data/red_propagada.txt --output data/string_interactions.tsv

# Análisis de la red con R
echo "Analizando propiedades de la red con R..."
Rscript propiedades_red.R data/string_interactions.tsv results/network_analysis_results.txt > logs/propiedades_red.log 2>&1

# Mover imágenes generadas por propiedades_red.R
if [[ -f "network_clustering_visualization.png" ]]; then
    echo "Guardando visualización de la red en 'images/'..."
    mv network_clustering_visualization.png images/network_clustering_visualization.png
fi

# Verificar y mover el archivo de clusters de genes
if [[ -f "clusters_genes.txt" ]]; then
    echo "Guardando clusters de genes en 'results/'..."
    mv clusters_genes.txt results/clusters_genes.txt
else
    echo "Error: No se generó el archivo 'clusters_genes.txt'."
    exit 1
fi

# Análisis de enriquecimiento funcional (descomentar cuando esté listo)
# echo "Realizando análisis de enriquecimiento funcional..."
# python analisis_enriquecimiento.py results/clusters_genes.txt results/enrichment_results.txt

echo "Pipeline completado. Los resultados están en la carpeta 'results/'."
