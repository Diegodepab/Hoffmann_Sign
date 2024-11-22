import pandas as pd
from tqdm import tqdm
import requests
import argparse

parser = argparse.ArgumentParser(description='Mapeo de genes a proteínas String.')
parser.add_argument('--genes_iniciales', type=str, help='Archivo con los genes a mapear.')
parser.add_argument('--genes_mapeados', type=str, help='Archivo de salida con los genes mapeados.')
args = parser.parse_args()


# Cargar el archivo con los genes 
df = pd.read_csv(args.genes_iniciales, sep='\t')


# Crear un cliente de MyGene
def map_to_protein(gene_symbols, species='9606'):
    """Mapea un lote de IDs a nombres de proteínas usando la API de STRING."""
    string_api_url = "https://string-db.org/api/tsv/get_string_ids"
    headers = {'User-Agent': 'Python-script'}
    preferred_names = []
    for gene_symbol in gene_symbols:
        params = {
            'identifiers': gene_symbol,
            'species': species
        }
        try:
            response = requests.get(string_api_url, params=params, headers=headers)
            if response.ok:
                lines = response.text.strip().split('\n')
                for line in lines[1:]:
                    fields = line.split('\t')
                    preferred_names.append([fields[1]])# Mapeo de ID a nombre preferido
            else:
                print(f"Error en la consulta para {gene_symbol}: {response.text}")
        except requests.RequestException as e:
            print(f"Error de conexión para {gene_symbol}: {e}")
    return preferred_names



mapped_genes = map_to_protein(df[' name'])
print(f"Se mapearon {len(mapped_genes)} genes a nombres de proteínas.")
# flat list para manejabilidad
mapped_genes = [item for sublist in mapped_genes for item in sublist]
#quitamos el prefijo de los genes
mapped_genes = [gene.split(".")[1] for gene in mapped_genes]

# Unir los genes mapeados en una cadena
gene_str = " ".join(mapped_genes)

# Guardar los resultados en un archivo
with open(args.genes_mapeados, 'w') as f:
    f.write( gene_str)


print("El mapeo de genes se completó con éxito.")