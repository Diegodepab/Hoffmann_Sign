import pandas as pd
import requests

# Cargar datos de genes
# Supón que genes es un DataFrame o lista con varias listas de genes por cluster
# Ejemplo de entrada: [{'cluster1': ['gen1', 'gen2', ...]}, {'cluster2': ['gen3', 'gen4', ...]}, ...]

def enrichment_analysis_enrichr_keeg_func(genes_by_cluster, combined_score_threshold=40):
    # Enrichr API URLs
    add_list_url = "https://maayanlab.cloud/Enrichr/addList"
    enrich_url = "https://maayanlab.cloud/Enrichr/enrich"

    all_results = {}

    # Para cada cluster de genes
    for cluster_name, gene_symbols in genes_by_cluster.items():
        print(f"Análisis de enriquecimiento para el cluster: {cluster_name}")

        # Convierte los nombres de los genes a mayúsculas
        gene_symbols_upper = [gene.upper() for gene in gene_symbols]

        # Agregar genes a Enrichr
        gene_str = "\n".join(gene_symbols_upper)
        payload = {
            'list': (None, gene_str),
            'description': (None, None)
        }
        response = requests.post(add_list_url, files=payload)
        if response.status_code == 200:
            result = response.json()
            user_list_id = result['userListId']
            print(f"Genes de {cluster_name} agregados correctamente a Enrichr. ID de lista de usuario: {user_list_id}")
        else:
            print(f"Error al agregar los genes de {cluster_name} a Enrichr: {response.status_code}")
            continue

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
                print(f"\nResultados del análisis de enriquecimiento para {cluster_name}:")
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
                print(f"La biblioteca de conjuntos de genes esperada no se encontró en los resultados de enriquecimiento para {cluster_name}.")
        else:
            print(f"Error al recuperar los resultados de enriquecimiento para {cluster_name}: {response.status_code}")
    
    # Comparar los resultados entre clusters
    compare_clusters(all_results)

def compare_clusters(all_results):
    """
    Compara los resultados de enriquecimiento entre diferentes clusters de genes.
    """
    for cluster_name, results in all_results.items():
        print(f"\nResultados de enriquecimiento para el cluster {cluster_name}:")
        for result in results:
            term_name, pvalue, combined_score, involved_genes, num_genes = result
            print(f"Term: {term_name}, P-value: {pvalue}, Combined Score: {combined_score}, Genes: {involved_genes}, Number of Genes: {num_genes}")

# Ejemplo de cómo pasar las listas de genes por cluster
genes_by_cluster = {
    'cluster1': ['gen1', 'gen2', 'gen3'],
    'cluster2': ['gen4', 'gen5', 'gen6']
}

print("Análisis de enriquecimiento usando Enrichr y rutas KEGG")
enrichment_analysis_enrichr_keeg_func(genes_by_cluster)
