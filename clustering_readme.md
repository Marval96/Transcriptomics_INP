# **Clustering espacial**

Se realizó el clustering espacial de los datos con base en tres radios: 50, 100 y 150 px. Para cada uno de ellos se corrió un análisis de componentes principales sobre la matriz espacial, con el objetivo de observar el comportamiento de los datos.

* Los gráficos se encuentran en:

  * Varianza del PCA; PCA_variance_spatial_count_$RADIO.png
  * PC1 vs PC2; PCA_PC1_vs_PC2_$RADIO.png

* Las matrices u otros archivos de texto se encuentran en:

  * Varianza explicada por cada PC; pca_variance_explained_$RADIO.csv
  * Selección del número óptimo de PC; pca_pc_selection_criteria_$RADIO.csv
  * Matriz de cargas; pca_loadings_$RADIO.csv
  * Variables principales por PC; pca_top_loadings_$RADIO.csv

Posteriormente, se evaluó el valor *k* de K-means utilizando distintos valores, con el objetivo de determinar el punto en el que se obtiene el mejor clustering con base en las métricas de Silhouette, Davies–Bouldin y Calinski–Harabasz. Estas métricas se normalizaron y promediaron para obtener un puntaje combinado y así seleccionar un valor *k* capaz de reflejar un clustering adecuado. Para este análisis se empleó la matriz espacial con todos los PC incluidos (17 componentes principales).

* Los gráficos se encuentran en:

  * Evaluación del clustering; kmeans_elbow_analysis_pcs17_$RADIO.png

* Las matrices u otros archivos de texto se encuentran en:

  * Evaluación del clustering; kmeans_optimization_pcs17_$RADIO.csv

El siguiente paso fue realizar el clustering mediante K-means, utilizando el valor *k* obtenido en el paso anterior.

* Los gráficos se encuentran en:

  * Visualización espacial de los clústeres; spatial_scatter_kmeans_spatial_count_$RADIO_17pcs.png
  * Visualización espacial de los clústeres ordenados; spatial_scatter_ordered_pcs17_$RADIO.png
  * Composición de RCN por fenotipo; cluster_phenotype_heatmap_pcs17_$RADIO.png
  * Enriquecimiento de los fenotipos en cada RCN; cluster_phenotype_enrichment_heatmap_pcs17_$RADIO.png
  * Proporción de fenotipos en RCN; stacked_barplot_kmeans_spatial_count_pcs17_$RADIO.png

* Las matrices u otros archivos de texto se encuentran en:

  * Composición de RCN por fenotipo; cluster_phenotype_compositionpcs17_$RADIO.csv
  * Enriquecimiento de los fenotipos en cada RCN; cluster_phenotype_enrichment_pcs17_$RADIO.csv

Finalmente, se evaluó mediante un PCA el comportamiento de los RCN con base en su cantidad (números absolutos) de células por fenotipo y su proporción de células positivas a marcadores funcionales. Para estudiar la similitud o redundancia entre estos RCN, se realizó un agrupamiento jerárquico considerando todas las dimensiones de los datos y un punto de corte en el percentil 10 o cuartil 25, con el fin de determinar redundancia y su posible agrupamiento supervisado en un metaclúster.

* Los gráficos se encuentran en:

  * PCA; cluster_PCA_similarity_rad$RADIO.png
  * Varianza del PCA; PCA_variance_cluster_characterization_rad$RADIO.png
  * Agrupamiento jerárquico; cluster_dendrogram_rad$RADIO_sim$CUTOFF.png
  * PCA y dendrograma; combined_pca_dendrogram_rad$RADIO_sim$CUTOFF.png

* Las matrices u otros archivos de texto se encuentran en:

  * Matriz de caracterización; cluster_characterization_matrix_rad$RADIO.csv
  * Varianza explicada por cada PC; pca_variance_table_rad$RADIO.csv
  * Selección del número óptimo de PC; pca_pc_selection_criteria_clusters_rad$RADIO.csv
  * Matriz de cargas; pca_loadings_clusters_rad$RADIO.csv
  * Variables principales por PC; pca_top_loadings_clusters_rad$150.csv
  * Distancias entre clústeres similares; similar_cluster_pairs_rad$RADIO_sim$CUTOFF.csv
  * Distancia entre clústeres; cluster_distance_matrix_rad$RADIO_sim$CUTOFF.csv
