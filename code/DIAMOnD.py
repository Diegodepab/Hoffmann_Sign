# Basado en el código de Susan Dina Ghiassian, Joerg Menche & Albert-Laszlo Barabasi


import time
import networkx as nx
import numpy as np
import copy
import scipy.stats
from collections import defaultdict
import csv
import sys


# =============================================================================
def print_usage():

    print(' ')
    print('        usage: python3 DIAMOnD.py network_file seed_file n alpha(optional) outfile_name (optional)')
    print('        -----------------------------------------------------------------')
    print('        network_file : The edgelist must be provided as any delimiter-separated')
    print('                       table. Make sure the delimiter does not exit in gene IDs')
    print('                       and is consistent across the file.')
    print('                       The first two columns of the table will be')
    print('                       interpreted as an interaction gene1 <==> gene2')
    print('        seed_file    : table containing the seed genes (if table contains')
    print('                       more than one column they must be tab-separated;')
    print('                       the first column will be used only)')
    print('        n            : desired number of DIAMOnD genes, 200 is a reasonable')
    print('                       starting point.')
    print('        alpha        : an integer representing weight of the seeds,default')
    print('                       value is set to 1')
    print('        outfile_name : results will be saved under this file name')
    print('                       by default the outfile_name is set to "first_n_added_nodes_weight_alpha.txt"')
    print(' ')


def check_input_style(input_list):
    try:
        network_edgelist_file = input_list[1]
        seeds_file = input_list[2]
        max_number_of_added_nodes = int(input_list[3])
    # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return

    alpha = 1
    outfile_name = 'first_%d_added_nodes_weight_%d.txt' % (max_number_of_added_nodes, alpha)

    if len(input_list) == 5:
        try:
            alpha = int(input_list[4])
            outfile_name = 'first_%d_added_weight_%d.txt' % (max_number_of_added_nodes, alpha)
        except:
            outfile_name = input_list[4]

    if len(input_list) == 6:
        try:
            alpha = int(input_list[4])
            outfile_name = input_list[5]
        except:
            print_usage()
            sys.exit(0)
            return
    return network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name


def read_input(network_file, seed_file):
    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file, 'r'):
        if line[0] == '#':
            continue
        else:
            dialect = sniffer.sniff(line)
            line_delimiter = dialect.delimiter
            break
    if line_delimiter == None:
        print
        'network_file format not correct'
        sys.exit(0)

    G = nx.Graph()
    for line in open(network_file, 'r'):
        if line[0] == '#':
            continue
        line_data = line.strip().split(line_delimiter)
        node1 = line_data[0]
        node2 = line_data[1]
        G.add_edge(node1, node2)

    seed_genes = set()
    for line in open(seed_file, 'r'):
        if line[0] == '#':
            continue
        line_data = line.strip().split('\t')
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)

    return G, seed_genes

def compute_all_gamma_ln(N):
    gamma_ln = {}
    for i in range(1, N + 1):
        gamma_ln[i] = scipy.special.gammaln(i)

    return gamma_ln

def logchoose(n, k, gamma_ln):
    if n - k + 1 <= 0:
        return scipy.infty
    lgn1 = gamma_ln[n + 1]
    lgk1 = gamma_ln[k + 1]
    lgnk1 = gamma_ln[n - k + 1]
    return lgn1 - [lgnk1 + lgk1]

def gauss_hypergeom(x, r, b, n, gamma_ln):
    return np.exp(logchoose(r, x, gamma_ln) +
                  logchoose(b, n - x, gamma_ln) -
                  logchoose(r + b, n, gamma_ln))

def pvalue(kb, k, N, s, gamma_ln):
    p = 0.0
    for n in range(kb, k + 1):
        if n > s:
            break
        prob = gauss_hypergeom(n, s, N - s, k, gamma_ln)
        # print prob
        p += prob

    if p > 1:
        return 1
    else:
        return p



def get_neighbors_and_degrees(G):
    neighbors, all_degrees = {}, {}
    for node in G.nodes():
        nn = set(G.neighbors(node))
        neighbors[node] = nn
        all_degrees[node] = G.degree(node)

    return neighbors, all_degrees

def reduce_not_in_cluster_nodes(all_degrees, neighbors, G, not_in_cluster, cluster_nodes, alpha):
    reduced_not_in_cluster = {}
    kb2k = defaultdict(dict)
    for node in not_in_cluster:

        k = all_degrees[node]
        kb = 0
        for neighbor in neighbors[node]:
            if neighbor in cluster_nodes:
                kb += 1
        k += (alpha - 1) * kb
        kb += (alpha - 1) * kb
        kb2k[kb][k] = node
    k2kb = defaultdict(dict)
    for kb, k2node in kb2k.items():
        min_k = min(k2node.keys())
        node = k2node[min_k]
        k2kb[min_k][kb] = node

    for k, kb2node in k2kb.items():
        max_kb = max(kb2node.keys())
        node = kb2node[max_kb]
        reduced_not_in_cluster[node] = (max_kb, k)

    return reduced_not_in_cluster


# ======================================================================================
#   C O R E    A L G O R I T H M
# ======================================================================================
def diamond_iteration_of_first_X_nodes(G, S, X, alpha):
    """
    Parameters:
    ----------
    - G:     graph
    - S:     seeds
    - X:     the number of iterations, i.e only the first X gened will be
             pulled in
    - alpha: seeds weight
    Returns:
    --------

    - added_nodes: ordered list of nodes in the order by which they
      are agglomerated. Each entry has 4 info:
      * name : dito
      * k    : degree of the node
      * kb   : number of +1 neighbors
      * p    : p-value at agglomeration
    """

    N = G.number_of_nodes()

    added_nodes = []

    # ------------------------------------------------------------------
    # Setting up dictionaries with all neighbor lists
    # and all degrees
    # ------------------------------------------------------------------
    neighbors, all_degrees = get_neighbors_and_degrees(G)

    # ------------------------------------------------------------------
    # Setting up initial set of nodes in cluster
    # ------------------------------------------------------------------

    cluster_nodes = set(S)
    not_in_cluster = set()
    s0 = len(cluster_nodes)

    s0 += (alpha - 1) * s0
    N += (alpha - 1) * s0

    # ------------------------------------------------------------------
    # precompute the logarithmic gamma functions
    # ------------------------------------------------------------------
    gamma_ln = compute_all_gamma_ln(N + 1)

    # ------------------------------------------------------------------
    # Setting initial set of nodes not in cluster
    # ------------------------------------------------------------------
    for node in cluster_nodes:
        not_in_cluster |= neighbors[node]
    not_in_cluster -= cluster_nodes

    # ------------------------------------------------------------------
    #
    # M A I N     L O O P
    #
    # ------------------------------------------------------------------

    all_p = {}

    while len(added_nodes) < X:

        # ------------------------------------------------------------------
        #
        # Going through all nodes that are not in the cluster yet and
        # record k, kb and p
        #
        # ------------------------------------------------------------------

        info = {}

        pmin = 10
        next_node = 'nix'
        reduced_not_in_cluster = reduce_not_in_cluster_nodes(all_degrees,
                                                             neighbors, G,
                                                             not_in_cluster,
                                                             cluster_nodes, alpha)

        for node, kbk in reduced_not_in_cluster.items():
            # Getting the p-value of this kb,k
            # combination and save it in all_p, so computing it only once!
            kb, k = kbk
            try:
                p = all_p[(k, kb, s0)]
            except KeyError:
                p = pvalue(kb, k, N, s0, gamma_ln)
                all_p[(k, kb, s0)] = p

            # recording the node with smallest p-value
            if p < pmin:
                pmin = p
                next_node = node

            info[node] = (k, kb, p)

        # ---------------------------------------------------------------------
        # Adding node with smallest p-value to the list of aaglomerated nodes
        # ---------------------------------------------------------------------
        added_nodes.append((next_node,
                            info[next_node][0],
                            info[next_node][1],
                            info[next_node][2]))

        # Updating the list of cluster nodes and s0
        cluster_nodes.add(next_node)
        s0 = len(cluster_nodes)
        not_in_cluster |= (neighbors[next_node] - cluster_nodes)
        not_in_cluster.remove(next_node)

    return added_nodes


# ===========================================================================
#
#   M A I N    D I A M O n D    A L G O R I T H M
#
# ===========================================================================
def DIAMOnD(G_original, seed_genes, max_number_of_added_nodes, alpha, outfile=None):
    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    if len(disease_genes) != len(seed_genes):
        print("DIAMOnD(): ignoring %s of %s seed genes that are not in the network" % (
            len(seed_genes - all_genes_in_network), len(seed_genes)))

    # 2. agglomeration algorithm.
    added_nodes = diamond_iteration_of_first_X_nodes(G_original,
                                                     disease_genes,
                                                     max_number_of_added_nodes, alpha)
    # 3. saving the results
    with open(outfile, 'w') as fout:

        fout.write('\t'.join(['#rank', 'DIAMOnD_node', 'p_hyper']) + '\n')
        rank = 0
        for DIAMOnD_node_info in added_nodes:
            rank += 1
            DIAMOnD_node = DIAMOnD_node_info[0]
            p = float(DIAMOnD_node_info[3])

            fout.write('\t'.join(map(str, ([rank, DIAMOnD_node, p]))) + '\n')

    return added_nodes

if __name__ == '__main__':
    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] file providing the network in the form of an edgelist
    #     (tab-separated table, columns 1 & 2 will be used)
    #
    # [2] file with the seed genes (if table contains more than one
    #     column they must be tab-separated; the first column will be
    #     used only)
    #
    # [3] number of desired iterations
    #
    # [4] (optional) seeds weight (integer), default value is 1
    # [5] (optional) name for the results file

    # check if input style is correct
    input_list = sys.argv
    network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name = check_input_style(input_list)

    # read the network and the seed genes:
    G_original, seed_genes = read_input(network_edgelist_file, seeds_file)

    # run DIAMOnD
    added_nodes = DIAMOnD(G_original,
                          seed_genes,
                          max_number_of_added_nodes, alpha,
                          outfile=outfile_name)

    print("\n results have been saved to '%s' \n" % outfile_name)
