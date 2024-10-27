import pandas as pd
import requests

# Cargamos los genes en un DataFrame
genes = pd.read_csv("genes.tsv", sep="\t", header=None, names=["id", "name"])

def enrichment_analysis_enrichr_keeg_func(gene_symbols, combined_score_threshold=40):
    # Define la lista de genes y conviértelos a mayúsculas
    gene_symbols_upper = [gene.upper() for gene in gene_symbols]

    # Enrichr API URLs
    add_list_url = "https://maayanlab.cloud/Enrichr/addList"
    enrich_url = "https://maayanlab.cloud/Enrichr/enrich"

    # Agregar genes a Enrichr
    gene_str = "\n".join(gene_symbols_upper)  # Enrichr requiere símbolos de genes separados por nueva línea

    payload = {
        'list': (None, gene_str),
        'description': (None, None)
    }
    response = requests.post(add_list_url, files=payload)
    if response.status_code == 200:
        result = response.json()
        user_list_id = result['userListId']
        print(f"Successfully added genes to Enrichr. User List ID: {user_list_id}")
    else:
        print(f"Error adding genes to Enrichr: {response.status_code}")
        return

    # Realizar análisis de enriquecimiento
    gene_set_library = 'KEGG_2019_Human'  # Ejemplo: rutas KEGG
    params = {
        'userListId': user_list_id,
        'backgroundType': gene_set_library
    }

    response = requests.get(enrich_url, params=params)
    if response.status_code == 200:
        enrichment_results = response.json()
        if gene_set_library in enrichment_results:
            print("\nEnrichment Analysis Results (Top Results):")
            filtered_results = []
            for result in enrichment_results[gene_set_library]:
                term_name = result[1]  # Nombre del término
                pvalue = result[2]     # P-valor
                combined_score = result[4]  # Puntuación combinada
                involved_genes = result[5]  # Genes involucrados

                # Filtrar por combined score
                if combined_score > combined_score_threshold:
                    num_genes = len(involved_genes)  # Contar el número de genes involucrados directamente
                    filtered_results.append((term_name, pvalue, combined_score, involved_genes, num_genes))

            # Ordenar resultados por número de genes involucrados y luego por combined score (en orden descendente)
            filtered_results.sort(key=lambda x: (x[4], x[2]), reverse=True)

            with open('enrichment_results_Enrichr.txt', 'w') as file:  # Abrir un archivo para escribir resultados
                for result in filtered_results:
                    term_name, pvalue, combined_score, involved_genes, num_genes = result
                    output = f"Term: {term_name}, P-value: {pvalue}, Combined Score: {combined_score}, Genes: {involved_genes}, Number of Genes: {num_genes}\n"
                    print(output.strip())
                    file.writelines(output)  # Escribir el resultado en el archivo
                print(f"Análisis guardado en enrichment_results_Enrichr.txt")
        else:
            print("La biblioteca de conjuntos de genes esperada no se encontró en los resultados de enriquecimiento.")
    else:
        print(f"Error retrieving enrichment results: {response.status_code}")

print("Análisis de enriquecimiento usando Enrichr y rutas KEGG")
enrichment_analysis_enrichr_keeg_func(genes["name"])
