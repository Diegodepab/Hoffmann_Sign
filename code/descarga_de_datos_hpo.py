import requests

# URL base del HPO
HPO="HP:0031993"
HPO_GENES=f"https://ontology.jax.org/api/network/annotation/{HPO}/download/gene"
HPO_DISEASES=f"https://ontology.jax.org/api/network/annotation/{HPO}/download/disease"

def download_data(url, filename):
    response = requests.get(url)
    if response.status_code == 200:
        with open(filename, 'wb') as file:
            file.write(response.content)
        print(f"Data successfully downloaded to {filename}")
    else:
        print(f"Failed to download data. Status code: {response.status_code}")

# Descargar genes y enfermedades
download_data(HPO_GENES, "genes.tsv")
download_data(HPO_DISEASES, "diseases.tsv")