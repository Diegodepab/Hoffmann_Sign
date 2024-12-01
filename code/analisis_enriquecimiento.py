import os
import pandas as pd
import requests
import argparse

# Asegúrate de que los directorios necesarios existan
results_dir = '../results/enrichment_analysis'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

def enrichment_analysis_enrichr_func(genes_by_cluster, combined_score_threshold=40):
    # URLs de la API de Enrichr
    add_list_url = "https://maayanlab.cloud/Enrichr/addList"
    enrich_url = "https://maayanlab.cloud/Enrichr/enrich"

    # Bibliotecas de Enrichr para el análisis
    gene_set_libraries = [
        "KEGG_2019_Human",
        "Reactome_2021",
        "GO_Biological_Process_2021",
        "GO_Cellular_Component_2021",
        "GO_Molecular_Function_2021"
    ]

    all_results = {}

    # Para cada cluster de genes
    for cluster_name, gene_symbols in genes_by_cluster.items():
        gene_symbols_upper = [gene.upper() for gene in gene_symbols]
        gene_str = "\n".join(gene_symbols_upper)
        payload = {
            'list': (None, gene_str),
            'description': (None, f"Cluster {cluster_name}")
        }
        response = requests.post(add_list_url, files=payload)

        if response.status_code == 200:
            result = response.json()
            user_list_id = result.get('userListId', None)
            if user_list_id:
                for gene_set_library in gene_set_libraries:
                    params = {
                        'userListId': user_list_id,
                        'backgroundType': gene_set_library
                    }
                    response = requests.get(enrich_url, params=params)

                    if response.status_code == 200:
                        enrichment_results = response.json()
                        if gene_set_library in enrichment_results:
                            grouped_results = []
                            for result in enrichment_results[gene_set_library]:
                                term_name = result[1]
                                pvalue = result[2]
                                combined_score = result[4]
                                involved_genes = result[5]  # Lista de genes

                                # Filtrar por puntuación combinada y términos asociados a más de un gen
                                if combined_score > combined_score_threshold and len(involved_genes) > 1:
                                    grouped_results.append((term_name, pvalue, combined_score, ", ".join(involved_genes)))

                            # Ordenar los resultados por puntuación combinada (en orden descendente)
                            grouped_results.sort(key=lambda x: x[2], reverse=True)

                            # Limitar a 15 términos como máximo
                            if len(grouped_results) > 15:
                                grouped_results = grouped_results[:15]

                            # Guardar solo si hay términos válidos
                            if grouped_results:
                                all_results.setdefault(cluster_name, {})[gene_set_library] = grouped_results
        else:
            continue

    # Guardar los resultados agrupados solo si hay resultados válidos
    if all_results:
        save_results_by_category(all_results, results_dir)

def save_results_by_category(all_results, results_dir):
    """
    Guarda los resultados del análisis en archivos CSV agrupados por categoría.
    """
    for cluster_name, categories in all_results.items():
        for category, results in categories.items():
            # Crear un DataFrame para cada categoría
            df = pd.DataFrame(results, columns=["Term", "P-value", "Combined Score", "Genes"])
            # Guardar en la carpeta de resultados
            sanitized_category = category.replace(" ", "_").replace("/", "_")
            output_file = f"{results_dir}/{cluster_name}_{sanitized_category}_enrichment_results.csv"
            df.to_csv(output_file, index=False)

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
            cluster_name, genes_str = line.strip().split(":")
            genes = [gene.strip() for gene in genes_str.split(",")]
            genes_by_cluster[cluster_name] = genes
    return genes_by_cluster

if __name__ == "__main__":
    # Configuración de los argumentos de la línea de comandos
    parser = argparse.ArgumentParser(description="Análisis de enriquecimiento de genes usando Enrichr.")
    parser.add_argument('file', help="Ruta del archivo que contiene los genes por cluster.")
    args = parser.parse_args()

    # Cargar los genes desde el archivo especificado
    genes_by_cluster = load_genes_from_file(args.file)

    # Iniciar el análisis de enriquecimiento
    enrichment_analysis_enrichr_func(genes_by_cluster)
