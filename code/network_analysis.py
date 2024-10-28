# network_analysis.py
import igraph as ig

def analyze_network(interaction_file, output_file):
    # Leer interacciones de proteínas
    edges = []
    with open(interaction_file) as file:
        next(file)  # Saltar la cabecera
        for line in file:
            fields = line.strip().split("\t")
            protein1 = fields[2]
            protein2 = fields[3]
            score = float(fields[5])
            edges.append((protein1, protein2, score))

    # Crear grafo en iGraph
    g = ig.Graph.TupleList(edges, weights=True, directed=False)

    # Clustering con el método de Louvain
    clusters = g.community_multilevel(weights="weight")
    print(f"Total de clústeres encontrados: {len(clusters)}")

    # Análisis de métricas topológicas
    degree_centrality = g.degree()
    betweenness_centrality = g.betweenness()
    closeness_centrality = g.closeness()

    # Imprimir las métricas de los primeros nodos
    for vertex in range(10):
        print(f"Nodo {g.vs[vertex]['name']}: Degree={degree_centrality[vertex]}, "
              f"Betweenness={betweenness_centrality[vertex]}, Closeness={closeness_centrality[vertex]}")

    # Guardar la red y los clústeres
    with open(output_file, "w") as file:
        for cluster_id, cluster in enumerate(clusters):
            genes_in_cluster = [g.vs[vertex]["name"] for vertex in cluster]
            file.write(f"Clúster {cluster_id + 1}: {', '.join(genes_in_cluster)}\n")
        print(f"Clústeres guardados en '{output_file}'")

# Ejecutar la función si el archivo se ejecuta directamente
if __name__ == "__main__":
    analyze_network("string_interactions.tsv", "network_clusters.txt")
