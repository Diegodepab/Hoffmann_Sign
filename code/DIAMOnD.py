import networkx as nx
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from collections import defaultdict

# ------------------- FUNCIONES DE DIAMOND -------------------
def read_input(network_file, seed_file):
    """Lee la red y la lista de genes de semilla desde archivos externos, y considera el puntaje combinado como peso."""
    # Leer la red con pesos
    G = nx.Graph()
    with open(network_file, 'r') as f:
        for line in f:
            if line.startswith("#") or "combined_score" in line:
                continue  # Omitir líneas de encabezado o comentarios
            node1, node2, score = line.strip().split(',')
            score = int(score)  # Convertimos el score a entero
            G.add_edge(node1, node2, weight=score)

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
    max_index = len(gamma_ln) - 1
    if r + b > max_index or x > max_index or r < x or b < (n - x):
        return 0  # Retorno seguro para valores fuera de rango

    return np.exp(gamma_ln[r] - (gamma_ln[x] + gamma_ln[r - x]) +
                  gamma_ln[b] - (gamma_ln[n - x] + gamma_ln[b - n]) -
                  gamma_ln[r + b])

def pvalue(kb, k, N, s, gamma_ln):
    """Calcula el p-valor para un nodo con kb enlaces a semillas y k enlaces totales."""
    return sum(gauss_hypergeom(n, s, N - s, k, gamma_ln) for n in range(kb, k + 1))

def diamond_iteration_of_first_X_nodes(G, S, X, alpha=1):
    """Ejecuta una iteración del algoritmo DIAMOnD considerando el puntaje combinado."""
    added_nodes = []
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes}
    degrees = {node: G.degree(node) for node in G.nodes}
    cluster_nodes = set(S)
    gamma_ln = compute_all_gamma_ln(len(G.nodes))

    while len(added_nodes) < X:
        min_p = float('inf')
        next_node = None
        for node in set(G.nodes) - cluster_nodes:
            k = degrees[node]  # Grado total del nodo
            kb = sum(1 for neighbor in neighbors[node] if neighbor in cluster_nodes)  # Grado en semillas
            p = pvalue(kb, k, len(G.nodes), len(cluster_nodes), gamma_ln)
            weight_sum = sum(G[node][neighbor]['weight'] for neighbor in neighbors[node] if neighbor in cluster_nodes)

            if p < min_p or (p == min_p and weight_sum > sum(G[next_node][n]['weight'] for n in neighbors[next_node] if n in cluster_nodes)):
                min_p = p
                next_node = node

        if next_node:
            added_nodes.append(next_node)
            cluster_nodes.add(next_node)

    return added_nodes



# ------------------- EJECUCIÓN -------------------
network_file = "red_diamond.txt"  # Cambia por la ruta de tu archivo
seed_file = "genes.txt"    # Cambia por la ruta de tu archivo
num_diamond_genes = 35          # Número de genes a agregar
alpha = 1

# Leer red y genes de semilla
G, seed_genes = read_input(network_file, seed_file)

# Ejecutar DIAMOnD
diamond_genes = diamond_iteration_of_first_X_nodes(G, seed_genes, num_diamond_genes, alpha)
print(f"Genes de semilla: {seed_genes}")
print(f"Genes de DIAMOnD: {diamond_genes}")

