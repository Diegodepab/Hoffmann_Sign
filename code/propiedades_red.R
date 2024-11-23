# Instalar argparse si no está instalado
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

# Verificar e instalar librerías necesarias
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)

# Configurar argparse para manejar los argumentos
parser <- ArgumentParser(description = "Análisis de propiedades de red")
parser$add_argument("input_file", help = "Ruta al archivo de la red propagada")
parser$add_argument("output_folder", help = "Ruta al folder donde se guardarán los resultados")
args <- parser$parse_args()

input_file <- args$input_file
output_folder <- args$output_folder

# Verificar si el archivo de entrada existe
if (!file.exists(input_file)) {
  stop("El archivo de entrada especificado no existe: ", input_file)
}

# Verificar o crear el folder de salida
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Cargar la red de interacciones
interaction_network <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convertir a un objeto igraph (grafo no dirigido)
g <- graph_from_data_frame(d = interaction_network[, c("preferredName_A", "preferredName_B")], directed = FALSE)

# Archivo de resumen
summary_file <- file.path(output_folder, "network_summary.txt")
sink(summary_file)
cat("===== ANÁLISIS DE LA RED =====\n\n")

# 1. Número de nodos y aristas
num_nodes <- gorder(g)
num_edges <- gsize(g)
cat("1. Número de nodos y aristas:\n")
cat("Número de nodos:", num_nodes, "\n")
cat("Número de aristas:", num_edges, "\n\n")

# 2. Grado promedio
degree_values <- degree(g)
average_degree <- mean(degree_values)
cat("2. Grado promedio:\n")
cat("Grado promedio de todos los nodos:", average_degree, "\n\n")

# 3. Centralidad de cercanía
centrality_closeness <- closeness(g)
top_closeness <- sort(centrality_closeness, decreasing = TRUE)[1:5]
cat("3. Centralidad de cercanía:\n")
cat("Top 5 nodos con mayor centralidad de cercanía:\n")
for (nodo in names(top_closeness)) {
  cat(" -", nodo, ":", round(top_closeness[nodo], 4), "\n")
}
cat("\n")

# 4. Centralidad de intermediación
betweenness_centrality <- betweenness(g)
top_betweenness <- sort(betweenness_centrality, decreasing = TRUE)[1:5]
cat("4. Centralidad de intermediación:\n")
cat("Top 5 nodos con mayor centralidad de intermediación:\n")
for (nodo in names(top_betweenness)) {
  cat(" -", nodo, ":", round(top_betweenness[nodo], 4), "\n")
}
cat("\n")

# 5. Modularidad
comunidades <- cluster_louvain(g)
modularity_value <- modularity(comunidades)
sizes_comunidades <- sizes(comunidades)
cat("5. Modularidad de la red:\n")
cat("Valor de modularidad:", modularity_value, "\n")
cat("Tamaño de las comunidades:\n")
print(sizes_comunidades)
cat("\n")

# 6. Densidad de la red
sparse_density <- edge_density(g)
cat("6. Densidad de la red:\n")
cat("Densidad de la red:", sparse_density, "\n\n")

# Cerrar archivo de resumen
sink()

# 7. Coeficiente de clustering local y visualización
vertex_clustering <- transitivity(g, type = "local")
V(g)$color <- ifelse(vertex_clustering > 0, "red", "lightblue")

visualization_file <- file.path(output_folder, "network_clustering_visualization.png")
png(visualization_file, width = 800, height = 800)
set.seed(321)
plot(g, vertex.size = 5, vertex.label = NA, main = "Coeficiente de Clustering Local")
dev.off()

# Guardar genes por cluster
genes_por_cluster <- split(V(g)$name, membership(comunidades))
clusters_file <- file.path(output_folder, "clusters_genes.txt")
file_conn <- file(clusters_file, "w")
for (cluster_id in names(genes_por_cluster)) {
  cat("Cluster", cluster_id, ":", paste(genes_por_cluster[[cluster_id]], collapse = ", "), "\n", file = file_conn)
}
close(file_conn)

# Mensaje de finalización
cat("Todos los resultados se han guardado en el folder:", output_folder, "\n")
