# Cargar librerías necesarias
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)

# Ruta al archivo de red
file_path <- "string_network_filtered_hugo-400.tsv"
output_file <- "network_analysis_results.txt"  # Archivo de salida

# Cargar la red de interacciones
interaction_network <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convertir a un objeto igraph (grafo no dirigido)
g <- graph_from_data_frame(d = interaction_network, directed = FALSE)

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
betweenness_normalized <- (betweenness_centrality - min(betweenness_centrality)) / 
  (max(betweenness_centrality) - min(betweenness_centrality))
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


# Calcular el coeficiente de clustering local por nodo en la red
# El coeficiente de clustering mide cuán cercanos están los vecinos de un nodo entre sí.
cat("7. Calculando el coeficiente de clustering local por nodo...\n")
vertex_clustering <- transitivity(g, type = "local")

# Asignar colores a los nodos según su coeficiente de clustering: rojo si el coeficiente es mayor que 0, azul claro si es 0
V(g)$color <- ifelse(vertex_clustering > 0, "red", "lightblue")

# Visualizar la red, donde los nodos con coeficiente de clustering positivo se muestran en rojo
# Esto resalta los nodos que están bien conectados dentro de su vecindario.
cat("8. Visualizando la red con el coeficiente de clustering local...\n")
set.seed(321)  # Fijar semilla para reproducibilidad
plot(g, vertex.size = 5, vertex.label = NA, main = "Coeficiente de Clustering Local")

# Obtener los clusters de la red utilizando el algoritmo de Louvain
# El algoritmo de Louvain busca dividir la red en comunidades densamente conectadas.
cat("9. Detectando clusters en la red utilizando el algoritmo de Louvain...\n")
clusters <- cluster_louvain(g)

# Obtener la membresía de cada nodo (es decir, a qué cluster pertenece cada nodo)
memberships <- membership(clusters)

# Asegurarnos de que los nombres de los nodos están correctamente definidos
nombres_nodos <- V(g)$name

# Definir la ruta del archivo de salida para guardar los resultados de los clusters
archivo_salida <- file.path("..", "results", "clusters_genes.txt")

# Abrir el archivo en modo de escritura
file_conn <- file(archivo_salida, "w")

# Organizar los genes por cada cluster en una lista estructurada
genes_por_cluster <- split(nombres_nodos, memberships)

# Escribir los genes de cada cluster en el archivo de salida
cat("10. Escribiendo los genes por cluster en el archivo de salida...\n")
for (cluster_id in names(genes_por_cluster)) {
  # Formato de salida: Cluster ID seguido por los genes en ese cluster
  cat("Cluster", cluster_id, ":", paste(genes_por_cluster[[cluster_id]], collapse = ", "), "\n", file = file_conn)
}

# Cerrar el archivo después de escribir
close(file_conn)

# Informar al usuario que la salida ha sido guardada correctamente
cat("La salida se ha guardado en:", archivo_salida, "\n")

# Realizar un análisis de enriquecimiento GO para cada cluster
# Esto nos ayudará a entender las funciones biológicas asociadas a los genes en cada cluster
cat("11. Realizando análisis de enriquecimiento GO para cada cluster...\n")
enrichment_results_clusters <- list()

# Realizar el análisis GO para cada grupo de genes en cada cluster
for (i in seq_along(genes_por_cluster)) {
  # Realizar enriquecimiento GO para los genes de cada cluster
  enrichment_results_clusters[[i]] <- enrichGO(
    gene          = genes_por_cluster[[i]],   # Genes del cluster actual
    OrgDb         = org.Hs.eg.db,              # Base de datos de anotación de genes humanos
    keyType       = "SYMBOL",                  # Tipo de identificador: símbolos de genes
    ont           = "BP",                      # Ontología de proceso biológico (BP)
    pAdjustMethod = "BH",                      # Método de ajuste de p-valor (Benjamini-Hochberg)
    pvalueCutoff  = 0.05,                      # Umbral de significancia para el p-valor
    qvalueCutoff  = 0.2                        # Umbral para el q-valor
  )
}

# Visualización de los resultados de enriquecimiento para un cluster seleccionado
cat("12. Visualizando los resultados de enriquecimiento GO para el cluster seleccionado...\n")
cluster_a_visualizar <- 1  # Definir el cluster a visualizar (por ejemplo, el primer cluster)

# Verificar que el cluster seleccionado sea válido y visualizar los resultados
if (cluster_a_visualizar %in% seq_along(enrichment_results_clusters)) {
  # Visualizar los resultados de enriquecimiento GO para el cluster seleccionado
  dotplot(enrichment_results_clusters[[cluster_a_visualizar]], 
          showCategory = 10,  # Mostrar las 10 categorías más significativas
          title = paste("GO Enrichment for Cluster", cluster_a_visualizar))
} else {
  cat("El número del cluster ingresado no es válido.\n")
}

# Guardar los resultados de enriquecimiento GO para el cluster seleccionado
cat("13. Guardando los resultados de enriquecimiento GO para el cluster seleccionado...\n")
cluster_a_guardar <- 1  # Definir el cluster a guardar (por ejemplo, el primer cluster)

# Verificar que el cluster sea válido antes de intentar guardar los resultados
if (cluster_a_guardar %in% seq_along(enrichment_results_clusters)) {
  # Definir la ruta para guardar los resultados en la carpeta 'results'
  archivo_salida <- file.path("..", "results", paste0("enrichment_cluster_", cluster_a_guardar, ".csv"))
  
  # Crear la carpeta 'results' si no existe
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Convertir los resultados de enriquecimiento a un data frame y guardarlo como archivo CSV
  write.csv2(
    as.data.frame(enrichment_results_clusters[[cluster_a_guardar]]),  # Convertir a data frame
    file = archivo_salida   # Guardar el archivo con el nombre adecuado
  )
} else {
  cat("El número del cluster ingresado no es válido.\n")
}

# Mensaje final indicando que el análisis ha concluido
cat("===== FIN DEL ANÁLISIS =====\n")

# Cerrar el archivo de salida
sink()
