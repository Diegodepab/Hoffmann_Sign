

# Obtener el argumento de entrada
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
parser$add_argument("output_file", help = "Ruta al archivo de salida")

# Parsear los argumentos
args <- parser$parse_args()

input_file <- args$input_file
output_file <- args$output_file

# Verificar si el archivo de entrada existe
if (!file.exists(input_file)) {
  stop("El archivo de entrada especificado no existe: ", input_file)
}

# Cargar la red de interacciones
interaction_network <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convertir a un objeto igraph (grafo no dirigido)
g <- graph_from_data_frame(d = interaction_network[, c("preferredName_A", "preferredName_B")], directed = FALSE)

# Crear el archivo de salida
sink(output_file)

cat("===== ANÁLISIS DE LA RED =====\n\n")

# Número de nodos y aristas
num_nodes <- gorder(g)
num_edges <- gsize(g)
cat("1. Número de nodos y aristas:\n")
cat("Número de nodos:", num_nodes, "\n")
cat("Número de aristas:", num_edges, "\n\n")

# Grado promedio
degree_values <- degree(g)
average_degree <- mean(degree_values)
cat("2. Grado promedio:\n")
cat("Grado promedio de todos los nodos:", average_degree, "\n\n")

# Centralidad de cercanía
centrality_closeness <- closeness(g)
top_closeness <- sort(centrality_closeness, decreasing = TRUE)[1:5]
cat("3. Centralidad de cercanía:\n")
cat("Top 5 nodos con mayor centralidad de cercanía:\n")
for (nodo in names(top_closeness)) {
  cat(" -", nodo, ":", round(top_closeness[nodo], 4), "\n")
}
cat("\n")

# Centralidad de intermediación
betweenness_centrality <- betweenness(g)
top_betweenness <- sort(betweenness_centrality, decreasing = TRUE)[1:5]
cat("4. Centralidad de intermediación:\n")
cat("Top 5 nodos con mayor centralidad de intermediación:\n")
for (nodo in names(top_betweenness)) {
  cat(" -", nodo, ":", round(top_betweenness[nodo], 4), "\n")
}
cat("\n")

# Modularidad
comunidades <- cluster_louvain(g)
modularity_value <- modularity(comunidades)
sizes_comunidades <- sizes(comunidades)
cat("5. Modularidad de la red:\n")
cat("Valor de modularidad:", modularity_value, "\n")
cat("Tamaño de las comunidades:\n")
print(sizes_comunidades)
cat("\n")

# Densidad de la red
sparse_density <- edge_density(g)
cat("6. Densidad de la red:\n")
cat("Densidad de la red:", sparse_density, "\n\n")

# Coeficiente de clustering local
vertex_clustering <- transitivity(g, type = "local")
V(g)$color <- ifelse(vertex_clustering > 0, "red", "lightblue")
cat("7. Coeficiente de clustering local por nodo calculado.\n\n")

# Visualizar la red con colores basados en coeficiente de clustering
png("network_clustering_visualization.png", width = 800, height = 800)
set.seed(321)
plot(g, vertex.size = 5, vertex.label = NA, main = "Coeficiente de Clustering Local")
dev.off()
cat("8. Visualización de la red guardada en 'network_clustering_visualization.png'.\n\n")

# Genes por cluster
genes_por_cluster <- split(V(g)$name, membership(comunidades))
output_cluster_file <- "clusters_genes.txt"
file_conn <- file(output_cluster_file, "w")

cat("9. Genes por cluster:\n")
for (cluster_id in names(genes_por_cluster)) {
  cat("Cluster", cluster_id, ":", paste(genes_por_cluster[[cluster_id]], collapse = ", "), "\n", file = file_conn)
}
close(file_conn)
cat("Genes por cluster guardados en:", output_cluster_file, "\n\n")

cat("\n===== FIN DEL ANÁLISIS =====\n")

# Cerrar el archivo de salida
sink()

