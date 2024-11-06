import networkx as nx
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from collections import defaultdict

# ------------------- FUNCIONES DE DIAMOND -------------------
def read_input(network_file, seed_file):
    """Lee la red y la lista de genes de semilla desde archivos externos."""
    # Leer la red
    G = nx.Graph()
    with open(network_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            node1, node2 = line.strip().split(',')
            G.add_edge(node1, node2)

    # Leer los genes de semilla
    seed_genes = set()
    with open(seed_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            seed_genes.add(line.strip().split(',')[0])

    return G, seed_genes

def compute_all_gamma_ln(N):
    """Precomputa todas las funciones gamma logarítmicas."""
    gamma_ln = {i: scipy.special.gammaln(i) for i in range(1, N + 1)}
    return gamma_ln

def gauss_hypergeom(x, r, b, n, gamma_ln):
    """Distribución hipergeométrica."""
    max_index = len(gamma_ln) - 1  # El último índice permitido
    if r + b > max_index or x > max_index or r < x or b < (n - x):
        # Si los índices están fuera de rango, manejar el error
        return 0  # Un p-valor predeterminado, puedes cambiarlo a 1 si es necesario

    return np.exp(gamma_ln[r] - (gamma_ln[x] + gamma_ln[r - x]) +
                  gamma_ln[b] - (gamma_ln[n - x] + gamma_ln[b - n]) -
                  gamma_ln[r + b])

def pvalue(kb, k, N, s, gamma_ln):
    """Calcula el p-valor para un nodo con kb enlaces a semillas y k enlaces totales."""
    return sum(gauss_hypergeom(n, s, N - s, k, gamma_ln) for n in range(kb, k + 1))

def diamond_iteration_of_first_X_nodes(G, S, X, alpha=1):
    """Ejecuta una iteración del algoritmo DIAMOnD para obtener genes adicionales."""
    added_nodes = []
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes}
    degrees = {node: G.degree(node) for node in G.nodes}
    cluster_nodes = set(S)
    gamma_ln = compute_all_gamma_ln(len(G.nodes))  # Precomputar gamma_ln para el número total de nodos

    while len(added_nodes) < X:
        min_p = float('inf')
        next_node = None
        for node in set(G.nodes) - cluster_nodes:
            k = degrees[node]  # Grado total del nodo
            kb = sum((1 for neighbor in neighbors[node] if neighbor in cluster_nodes))  # Grado en semillas
            # Debug: imprimir valores
            print(f"Node: {node}, k: {k}, kb: {kb}, Cluster Size: {len(cluster_nodes)}")
            try:
                p = pvalue(kb, k, len(G.nodes), len(cluster_nodes), gamma_ln)
            except ValueError as e:
                print(f"Error calculando pvalue para nodo {node}: {e}")
                continue  # Salir del ciclo si hay error

            if p < min_p:
                min_p = p
                next_node = node

        if next_node:
            added_nodes.append(next_node)
            cluster_nodes.add(next_node)
    return added_nodes

# ------------------- VISUALIZACIÓN -------------------
def graficar_red_enriquecida(G, seed_genes, diamond_genes):
    """Crea y muestra el gráfico de la red enriquecida."""
    pos = nx.spring_layout(G)
    plt.figure(figsize=(10, 10))

    # Nodos de semillas
    nx.draw_networkx_nodes(G, pos, nodelist=seed_genes, node_color='lightblue', node_size=500, label="Seed Genes")

    # Nodos DIAMOnD
    nx.draw_networkx_nodes(G, pos, nodelist=diamond_genes, node_color='orange', node_size=300, label="DIAMOnD Genes")

    # Enlaces
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v) for u, v in G.edges() if u in (seed_genes | set(diamond_genes)) and v in (seed_genes | set(diamond_genes))], alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8)

    plt.legend(loc="best")
    plt.title("Red de Genes Enriquecida usando DIAMOnD")
    plt.show()

# ------------------- EJECUCIÓN -------------------
network_file = "red_diamond.txt"  # Cambia por la ruta de tu archivo
seed_file = "seed_ejemplo.txt"    # Cambia por la ruta de tu archivo
num_diamond_genes = 10            # Número de genes a agregar
alpha = 1

# Leer red y genes de semilla
G, seed_genes = read_input(network_file, seed_file)

# Ejecutar DIAMOnD
diamond_genes = diamond_iteration_of_first_X_nodes(G, seed_genes, num_diamond_genes, alpha)

# Graficar red enriquecida
graficar_red_enriquecida(G, seed_genes, diamond_genes)