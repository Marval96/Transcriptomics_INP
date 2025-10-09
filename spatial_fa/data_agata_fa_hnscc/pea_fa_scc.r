#---- Datos Generales ----------------------------------------------------------

# Tema: Analisis de enriquecimiento de vias.
# Objetivo: realizar un Analisis de Enriquecimineto de Vias (PEA) mediante ORA.
# Programador: Raul Valderrama.
# Version de R: 4.5
# Fecha: Mayo, 2025.

# Un análisis de sobrerepresentación determinia un valor estadistico 
# de la presencia de un set definido de genes (GO, KEEG, etc) 
# en un set de genes de interes.
#---------------------------------------------------------------------------------

#---- Paqueterias --------------------------------------------------------------

# Instalar paqueterias:

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("biomaRt")
#install.packages("wordcloud")
#install.packages("extrafont")
#BiocManager::install("ReactomePA")

# Activar paqueterias:

# Para genererar el objeto PEA..
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReactomePA)

# Para los graficos de PEA.
library(ggplot2)
library(enrichplot)
library(ggupset)
library(wordcloud)
library(RColorBrewer)
library(dplyr)

# Topologia
library(pathview) 
library(ggnewscale)

# Para convertir los "gen name" en Id de ENSMBLE.
library(biomaRt)

#---- Importar Base de Datos ---------------------------------------------------

# Establecer directorio de trabajo
setwd("./")
getwd()
list.files()

# La base de datos (db) contiene todos los genes expresados de manera diferencial. 
# El archivo debe contener: nombre del gen, ID ENSEMBLE y el valor del foldchange.
# En caso de no tener el ID, este se puede obtener a traves de BIOMART.
# Leer el archivo CSV con los genes
data <- read.csv("FA_SCC_filtered_down_genes.csv")
head(data)
dim(data)


# Up & down
#upregulated_genes <- data[data$log2FoldChange >= 0.58, ]
downregulated_genes <- data[data$log2FoldChange <= -0.58, ]
# Multiplicar log2fc por -1 para que queden como positivos
downregulated_genes$log2FoldChange <- downregulated_genes$log2FoldChange * -1


# Verificar nombres de columna
#head(upregulated_genes)
#colnames(upregulated_genes)
head(downregulated_genes)
colnames(downregulated_genes)

# Dimensiones de la db
dim(data)
dim(upregulated_genes)
dim(downregulated_genes)

#---- ID ENSEMBLE --------------------------------------------------------------

# Convertir los nombres de genes a Entrez IDs
genes_entrez <- bitr(#upregulated_genes$gene_symbol,
                     downregulated_genes$gene_symbol,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

# Extraer los IDs únicos
gene_list <- unique(genes_entrez$ENTREZID)
head(gene_list)

#---- Enriquecimiento: ORA -----------------------------------------------------

pea_ora <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",               # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)         # convierte ENTREZ a SYMBOL

# # Con Reactome
# pea_ora <- enrichPathway(
#     gene = gene_list,
#     organism = "human",         # O "Homo sapiens"
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE             # Convierte ENTREZ a símbolos
# )

nrow(pea_ora@result)

# Salvar resultados 
write.csv(as.data.frame(pea_ora), 
          file = "pae_go_bp_fa_scc_down.csv", row.names = FALSE)

# Nube de terminos
# Extraer la tabla de resultados del objeto enrichResult
df_go <- pea_ora@result
# Separar el numerador del GeneRatio (genes enriquecidos por término)
#numerador <- read.table(text = df_go$GeneRatio, sep = "/")[1]
numeradores <- sapply(strsplit(df_go$GeneRatio, "/"), function(x) as.numeric(x[1]))

# Crear data frame con la frecuencia (número de genes por término) y nombre del término
wcdf <- data.frame(
  term = df_go$Description,
  #freq = numerador$V1
  freq = numeradores
)

# Generar la nube de palabras
wordcloud(words = wcdf$term,
            freq = wcdf$freq,
            scale = c(2, 0.01),
            colors = brewer.pal(8, "Dark2"),
            max.words = 200)

# Salvar nube en tiff
tiff("wordcloud_go_bp_fa_scc_down_ora.tiff", width = 3000, height = 3000, res = 600)
wordcloud(words = wcdf$term,
          freq = wcdf$freq,
          scale = c(2, 0.01),
          colors = brewer.pal(8, "Dark2"),
          max.words = 200)
dev.off()



#----  Redundancia -------------------------------------------------------------

# Simplificar términos redundantes
nrow(pea_ora@result)
pea_ora <- simplify(pea_ora,
                    cutoff = .5,
                    by = "p.adjust",  # nolint
                    select_fun = min, 
                    measure = "Wang")
nrow(pea_ora@result)

# Filtrar con reactome
# Agrupar pathways por número de genes compartidos
# overlap <- pea_ora@result %>%
#   group_by(ID) %>%
#   summarize(genes = strsplit(geneID, "/"),
#             n_genes = lengths(genes))
# 
# # puedes filtrar por n_genes para mostrar solo los más específicos
# especificos <- overlap %>% filter(n_genes < 50)
# ids_especificos <- especificos$ID
# pea_ora_filtrado <- pea_ora
# pea_ora_filtrado@result <- pea_ora@result %>% filter(ID %in% ids_especificos)
# 
# 
# nrow(pea_ora_filtrado@result)



# Salvar resultados 
write.csv(as.data.frame(pea_ora), 
          file = "fa_scc_down_ora_redu_go_bp_results.csv", row.names = FALSE)

#---- Visualizacón--------------------------------------------------------------


dotplot(pea_ora, showCategory = 20)
tiff("dotplot_go_macs_uivc4_ora.tiff", width = 3000, height = 4000, res = 300)
dotplot(pea_ora, showCategory = 20)
dev.off()

p <- barplot(pea_ora, showCategory = 20)
p
pp <- p + scale_fill_viridis_c()
pp

tiff("barplot_fa_scc_down_ora.tiff", width = 6000, height = 7500, res = 600)
pp
dev.off()

# Traslape de genes en distintos terminos
#upsetplot(pea_ora)

# Topologia
# map_plot <- pairwise_termsim(pea_ora)
# map_plot <- emapplot(map_plot, showCategory = 20)
# map_plot

# Visualizar red de términos enriquecidos y genes asociados
# categorySize puede ser "pvalue" o "geneNum"
genes <- downregulated_genes$log2FoldChange
names(genes) <- downregulated_genes$gene_symbol
genes <- na.omit(genes)

net <- cnetplot(pea_ora, showCategory = 10,
                categorySize = "pvalue", 
                foldChange = genes) + scale_color_gradient(low = "blue", high = "red")

net

tiff("cnetplot_fa_scc_redu_reactome_ora.tiff", width = 3000, height = 3000, res = 300)
net
dev.off()


#---- Referencias: -------------------------------------------------------------

# https://rpubs.com/pranali018/enrichmentAnalysis
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/