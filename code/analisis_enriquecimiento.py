import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import argparse

# Asegúrate de que los directorios necesarios existan
if not os.path.exists('data/enrichment_analysis'):
    os.makedirs('data/enrichment_analysis')
if not os.path.exists('../results/enrichment_analysis'):
    os.makedirs('../results/enrichment_analysis')

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
                        filtered_results = []
                        for result in enrichment_results[gene_set_library]:
                            term_name = result[1]
                            pvalue = result[2]
                            combined_score = result[4]
                            involved_genes = result[5]

                            # Filtrar por puntuación combinada
                            if combined_score > combined_score_threshold:
                                num_genes = len(involved_genes)
                                filtered_results.append((term_name, pvalue, combined_score, involved_genes, num_genes))

                        # Ordenar resultados por número de genes involucrados y luego por combined score
                        filtered_results.sort(key=lambda x: (x[4], x[2]), reverse=True)

                        # Guardar los resultados para este cluster
                        all_results[cluster_name] = filtered_results
        else:
            continue

    # Guardar los resultados de enriquecimiento para cada cluster en un CSV y generar gráficos
    if all_results:
        save_results_and_plot(all_results)

def save_results_and_plot(all_results):
    """
    Guarda los resultados del análisis en archivos CSV y genera gráficos.
    """
    for cluster_name, results in all_results.items():
        # Guardar los resultados de enriquecimiento en un archivo CSV en 'data/enrichment_analysis'
        df = pd.DataFrame(results, columns=["Term", "P-value", "Combined Score", "Genes", "Number of Genes"])
        df.to_csv(f"data/enrichment_analysis/enrichment_results_{cluster_name}.csv", index=False)

        # Generar un gráfico de barras (p-value vs Term)
        plt.figure(figsize=(10, 6))
        top_results = df.head(10)  # Mostrar los 10 mejores resultados
        sns.barplot(x='P-value', y='Term', data=top_results, palette='viridis', hue='Term', legend=False)
        plt.title(f'Análisis de enriquecimiento para {cluster_name}')
        plt.xlabel('P-value')
        plt.ylabel('Term')
        plt.tight_layout()

        # Guardar el gráfico en 'results/enrichment_analysis'
        plt.savefig(f"../results/enrichment_analysis/enrichment_barplot_{cluster_name}.png")
        plt.close()

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
