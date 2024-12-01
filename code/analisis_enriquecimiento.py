import os
import pandas as pd
import requests
import argparse

# Asegúrate de que los directorios necesarios existan
results_dir = '../results/enrichment_analysis'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

def enrichment_analysis_enrichr_keeg_func(genes_by_cluster, combined_score_threshold=40):
    # URLs de la API de Enrichr
    add_list_url = "https://maayanlab.cloud/Enrichr/addList"
    enrich_url = "https://maayanlab.cloud/Enrichr/enrich"

    all_results = {}

    # Para cada cluster de genes
    for cluster_name, gene_symbols in genes_by_cluster.items():
        gene_symbols_upper = [gene.upper() for gene in gene_symbols]
        gene_str = "\n".join(gene_symbols_upper)
        payload = {
            'list': (None, gene_str),
            'description': (None, None)
        }
        response = requests.post(add_list_url, files=payload)

        if response.status_code == 200:
            result = response.json()
            user_list_id = result.get('userListId', None)
            if user_list_id:
                # Realizar análisis de enriquecimiento
                gene_set_library = 'KEGG_2019_Human'
                params = {
                    'userListId': user_list_id,
                    'backgroundType': gene_set_library
                }

                response = requests.get(enrich_url, params=params)

                if response.status_code == 200:
                    enrichment_results = response.json()

                    if gene_set_library in enrichment_results:
                        grouped_results = {}
                        for result in enrichment_results[gene_set_library]:
                            term_name = result[1]
                            pvalue = result[2]
                            combined_score = result[4]
                            involved_genes = result[5]

                            # Filtrar por puntuación combinada
                            if combined_score > combined_score_threshold:
                                term_category = categorize_term(term_name)
                                if term_category not in grouped_results:
                                    grouped_results[term_category] = []
                                grouped_results[term_category].append((term_name, pvalue, combined_score, involved_genes))

                        # Guardar los resultados agrupados para este cluster
                        all_results[cluster_name] = grouped_results
        else:
            continue

    # Guardar los resultados de enriquecimiento agrupados en tablas por categoría
    if all_results:
        save_grouped_results(all_results, results_dir)

def categorize_term(term_name):
    """
    Categoriza un término en función de su contenido.
    """
    if "pathway" in term_name.lower():
        return "Pathway"
    elif "cell" in term_name.lower():
        return "Cellular Function"
    elif "regulation" in term_name.lower():
        return "Regulation"
    else:
        return "Other"

def save_grouped_results(all_results, results_dir):
    """
    Guarda los resultados del análisis en archivos CSV agrupados por categoría.
    """
    for cluster_name, grouped_results in all_results.items():
        for category, results in grouped_results.items():
            # Crear un DataFrame para cada categoría
            df = pd.DataFrame(results, columns=["Term", "P-value", "Combined Score", "Genes"])
            # Guardar en la carpeta de resultados
            df.to_csv(f"{results_dir}/{cluster_name}_{category}_enrichment_results.csv", index=False)

def load_genes_from_file(file_path):
    """
    Carga los datos de genes desde un archivo de texto.
    El archivo debe tener el siguiente formato:
    cluster_name1: gene1, gene2, gene3
    cluster_name2: gene4, gene5, gene6
    """
    genes_by_cluster = {}

    with open(file_path, 'r') as f:
        for line in f:
            # El archivo debe tener una línea por cada cluster, donde los genes están separados por comas
            cluster_name, genes_str = line.strip().split(":")
            genes = [gene.strip() for gene in genes_str.split(",")]
            genes_by_cluster[cluster_name] = genes

    return genes_by_cluster

if __name__ == "__main__":
    # Configuración de los argumentos de la línea de comandos
    parser = argparse.ArgumentParser(description="Análisis de enriquecimiento de genes usando Enrichr y KEGG.")
    parser.add_argument('file', help="Ruta del archivo que contiene los genes por cluster.")
    args = parser.parse_args()

    # Cargar los genes desde el archivo especificado
    genes_by_cluster = load_genes_from_file(args.file)

    # Iniciar el análisis de enriquecimiento
    enrichment_analysis_enrichr_keeg_func(genes_by_cluster)
