#!/bin/bash

export R_LIB: ./R_packages

export PYTHON_LIB=./py_packages

mkdir -p data

# Descarga de datos de genes de HPO
wget https://ontology.jax.org/api/network/annotation/HP:0031993/download/gene -o /data/genes.tsv
# Pasamos de genes a string id
# python genes2string.py /data/genes.tsv /data/genes_string.tsv

# Desacarga de la red de String
wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz -o /data/network.txt

# Propagación de genes
# python dyamond.py /data/string.tsv /data/network.txt /data/red_propagada.txt

# Análisis de red
# python analisis.py /data/red_propagada.txt /data/clusters.txt

# Análisis de enriquecimiento funcional
# python analisis_enriquecimiento.py /data/clusters.txt 