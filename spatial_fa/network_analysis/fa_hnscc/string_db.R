# ---- STRING: redes de Interaacción Proteína-Proteína ----
# Este script es para obtener las interacciones entre proteínas reportadas
# en la base de datos de STRING. La plataforma web de STRING soporta un máximo 
# de 2000 porteínas. Este script permite superar este problema. 


# ---- Librería ----

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("STRINGdb")

install.packages("BiocManager")
BiocManager::install("STRINGdb")
library(STRINGdb)

# ----  Directorio de trabajo ----
setwd("./")
getwd()
list.files()

# ---- Manejo de datos ----

# Especifica el organismo; 
# (9606 = Homo sapiens); versión de la base de datos y score mínimo
string_db <- STRINGdb$new(
  version = "12",           # Última versión estable
  species = 9606,           # Homo sapiens
  score_threshold = 700,    # Score de confianza (400 es medio; puedes usar 700 para alto)
  input_directory = ""
)

# Listado de genes
list.files()

genes <- read.table("FA_SCC_filtered_up_genes.csv", 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
colnames(genes)
head(genes)
# Todos los genes son diferencialmente expresados
all(genes$padj < 0.05)
all(genes$log2FoldChange > 0.58)

# Obtener ID de string
mapped_genes <- string_db$map(genes, "gene_symbol", removeUnmappedRows = TRUE)
colnames(mapped_genes)

# Verifcar genes
# Filtra solo los que tienen STRING_id válido
mapped_genes_valid <- mapped_genes[!is.na(mapped_genes$STRING_id), ]
string_db$plot_network(mapped_genes_valid$STRING_id)

# Filtra NA y vacíos
ids <- mapped_genes$STRING_id
ids <- ids[!is.na(ids) & ids != ""]
ids <- unique(ids)   # Elimina duplicados

length(mapped_genes_valid$STRING_id)
head(mapped_genes_valid$STRING_id)

# Esto abre la red en el navegador
string_db$plot_network(mapped_genes$STRING_id)

# Obtener interacciones entre tus genes
interactions <- string_db$get_interactions(mapped_genes$STRING_id)
head(interactions)
# Paso previo a limpiar direccionalidad: eliminar duplicados exactos
interacciones <- unique(interactions)

# Ordenar los nodos para que A-B y B-A sean tratados igual
interacciones$node1_clean <- pmin(interacciones$from, interacciones$to)
interacciones$node2_clean <- pmax(interacciones$from, interacciones$to)

# Eliminar duplicados (B-A si ya existe A-B)
interacciones_nodir <- interacciones[!duplicated(interacciones[c("node1_clean", "node2_clean")]), ]

# Guardar solo las columnas necesarias
interacciones_final <- interacciones_nodir[, c("node1_clean", "node2_clean", "combined_score")]
colnames(interacciones_final) <- c("from", "to", "combined_score")

# Exportar como TSV compatible con Cytoscape o Excel
#write.table(interacciones_final, "interacciones_no_dirigidas.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Usar mapping inverso de STRING ID a símbolo de gen
genes_mapeados <- mapped_genes[, c("gene_name", "STRING_id")]
colnames(genes_mapeados) <- c("symbol", "STRING_id")
head(genes_mapeados)

# Hacer el merge para tener los nombres legibles
interacciones_final_named <- merge(interacciones_final, genes_mapeados, by.x="from", by.y="STRING_id")
interacciones_final_named <- merge(interacciones_final_named, genes_mapeados, by.x="to", by.y="STRING_id")

# Seleccionar columnas y renombrar
interacciones_final_named <- interacciones_final_named[, c("symbol.x", "symbol.y", "combined_score")]
colnames(interacciones_final_named) <- c("from", "to", "combined_score")

# Exportar
head(interacciones_final_named)
colnames(interacciones_final_named) <- c("#node1", "node2", "combined_score")
head(interacciones_final_named)
write.table(interacciones_final_named, "string_interactions_short.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#...............................................................................

# Leer archivo completo (tarda unos segundos, es grande)
list.files()
#all_links <- read.delim("9606.protein.links.full.v12.0.txt")
all_links <- read.table(gzfile("9606.protein.links.full.v12.0.txt"), 
                        header = TRUE, sep = " ", stringsAsFactors = FALSE)
# 2. Verifica que tiene las columnas correctas
head(all_links)
class(all_links)

# Filtrar solo las interacciones entre tus genes mapeados
ids <- mapped_genes$STRING_id
filtered_links <- all_links[all_links$protein1 %in% ids & all_links$protein2 %in% ids, ]

# Opcional: eliminar duplicados dirigidos
filtered_links$node1 <- pmin(filtered_links$protein1, filtered_links$protein2)
filtered_links$node2 <- pmax(filtered_links$protein1, filtered_links$protein2)
filtered_links <- filtered_links[!duplicated(filtered_links[c("node1", "node2")]), ]
head(filtered_links)
# Exportar con todas las columnas
write.table(filtered_links, "interacciones_completas_STRING.tsv", sep = "\t", quote = FALSE, row.names = FALSE)










