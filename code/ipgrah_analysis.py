import pandas as pd
import requests
import shutil

# Cargar los genes desde el archivo genes.tsv
genes_df = pd.read_csv("genes.tsv", sep="\t")
gene_names = genes_df[" name"].tolist()

# Parámetros de la solicitud
base_url = "https://string-db.org/api/image/network"
species_id = 9606  # Homo sapiens
required_score = 700  # Puntaje mínimo de interacción
add_white_nodes = 10  # Nodos extra para completar la red
network_type = "functional"

# Preparar los identificadores de genes para la solicitud
identifiers = "%0d".join(gene_names)  # Formato URL de StringDB

# Construir la URL completa para la solicitud a STRING
url_request = f"{base_url}?identifiers={identifiers}&species={species_id}&required_score={required_score}&add_white_nodes={add_white_nodes}&network_type={network_type}"

# Realizar la solicitud para obtener la imagen
print("Requesting network image from STRING: " + url_request)
response = requests.get(url_request, stream=True)

# Guardar la imagen de la red si la solicitud es exitosa
if response.status_code == 200:
    with open('gene_network.png', 'wb') as out_file:
        shutil.copyfileobj(response.raw, out_file)
    print("La imagen de la red de interacciones ha sido descargada como 'gene_network.png'")
else:
    print(f"Error en la solicitud: {response.status_code}, Mensaje: {response.text}")


