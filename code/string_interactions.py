# stringdb_interaction.py
import pandas as pd
import requests

def download_interactions(genes_file, output_file):
    # Cargar los genes desde el archivo genes.tsv
    genes_df = pd.read_csv(genes_file, sep="\t")
    gene_names = genes_df[" name"].tolist()

    # Preparar la URL para StringDB
    species_id = 9606  # Homo sapiens
    string_api_url = "https://string-db.org/api/tsv/network"
    genes_str = "%0d".join(gene_names)  # Concatenar genes con formato URL

    # Construir la solicitud a StringDB
    params = {
        "identifiers": genes_str,  # Genes separados por %0d (formato StringDB)
        "species": species_id,
        "required_score": 400,  # Puntaje mínimo (ajustable)
        "network_type": "functional"
    }

    # Realizar la solicitud GET a StringDB
    response = requests.get(string_api_url, params=params)
    if response.status_code == 200:
        with open(output_file, "w") as file:
            file.write(response.text)
        print(f"Red de interacciones descargada exitosamente en '{output_file}'")
    else:
        print(f"Error al descargar datos de StringDB: {response.status_code}")

# Ejecutar la función si el archivo se ejecuta directamente
if __name__ == "__main__":
    download_interactions("genes.tsv", "string_interactions.tsv")
