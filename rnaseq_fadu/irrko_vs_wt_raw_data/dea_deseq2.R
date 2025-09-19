#...............................................................................

# ---- Data ----
# Expresión Diferencial DESeq2
# Script para hacer expresión diferencial usando DESeq2 a partir de cuentas
# generadas con STAR.
# R version: 4.5.0


#...............................................................................

# ---- Establercer el directorio de trabajo ----
# version
setwd("./")
getwd()
list.files()

#...............................................................................

# ---- Instalar librerías ----

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("biomaRt")
# BiocManager::install("apeglm")
# install.packages("readr")
#BiocManager::install('EnhancedVolcano')

library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(dplyr)
library(readr)
library(apeglm)
library(EnhancedVolcano)

#...............................................................................

# ---- Importar datos de STAR ----

# Import txt matrix from STAR
count_matrix <- read.table("matrix_counts_DESeq2.txt", header = TRUE, sep = "\t", 
                   row.names = 1)

# Verifica el formato
head(count_matrix)
colnames(count_matrix)
colnames(count_matrix) <- c("WT_1", "WT_2", "WT_3", 
                    "IrrKO_1", "IrrKO_2", "IrrKO_3")

head(count_matrix)

# Define el vector de condiciones
condition <- factor(c(rep("WT", 3), rep("IrrKO", 3)))

# Crea el data frame de diseño
colData <- data.frame(row.names = colnames(count_matrix), condition)
head(colData)

#...............................................................................

# ---- Crear el objeto DESeq ----
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition)

#...............................................................................

# ---- Filtrado de genes ----

# Filtrado de genes con múltiples condiciones

# Definir criterios de expresión mínima
min_counts <- 10
min_samples <- 3

# Calcular genes expresados en al menos `min_samples` muestras por cada condición
expressed_in_group <- lapply(levels(dds$condition), function(cond) {
  rowSums(counts(dds)[, dds$condition == cond] >= min_counts) >= min_samples
})

# Conservar genes que cumplen el criterio en al menos un grupo
keep <- Reduce("|", expressed_in_group)

# Aplicar el filtro
dds <- dds[keep, ]

#...............................................................................

# ---- Correr DESeq2 ----

dds <- DESeq(dds)

#...............................................................................

# ---- Datos Normalizados ----

# Obtener los conteos normalizados
normalized_counts <- counts(dds, normalized = TRUE)
# Convertir los conteos a un data frame para manipularlos
normalized_counts_df <- as.data.frame(normalized_counts)
# Guardar los conteos normalizados en un archivo CSV
write.csv(normalized_counts_df, file = "normalized_RNASEQ_FaDu.csv",
          row.names = TRUE)

# Aplicar transformación VST: mayor control de la varianza
#vsd <- varianceStabilizingTransformation(dds)
# Obtener los conteos transformados
#vsd_counts <- assay(vsd)
# Convertir a un data frame
#vsd_counts_df <- as.data.frame(vsd_counts)
# Guardar los conteos transformados en un archivo CSV
#write.csv(vsd_counts_df, file = "vsd_counts_all_samples_L1.csv", row.names = TRUE)

#...............................................................................

# ---- Análisis de Expresion diferencial ----
# confg ideal: results <- results(dds, contrast = c("condition", "treated", "control"))
# dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
res <- results(dds, contrast = c("condition", "IrrKO", "WT"), alpha = 0.05,
               lfcThreshold = 0.58)#, altHypothesis = "greaterAbs")
res
summary(res)

# Convertir los resultados a un data.frame
res_df <- as.data.frame(res)

# Añadir nombre de los genes
res_df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_df), 
                        keytype = "ENSEMBL", column = "SYMBOL")

# Remover valores NA
res_df <- na.omit(res_df)

# Identificar genes duplicados en los nombres de fila
duplicados <- duplicated(rownames(res_df))
# Mostrar si hay genes duplicados
any(duplicados)  # Devuelve TRUE si hay genes duplicados

# Guardar todos los genes (con y sin significancia)
write.csv(res_df, file = "all_genes_WT_vs_IrrKO.csv", row.names = TRUE)

# Filtrar los genes que tienen un valor de p.adj menor a 0.05
res_df<- res_df[res_df$padj < 0.05, ]

# Guardar los resultados como un archivo CSV
write.csv(res_df, file = "differential_expression_results_WT_vs_uIrrKO.csv",
          row.names = TRUE)
# Resumen
summary(res_df)

#...............................................................................

# ---- Volcano Plot ----

# Activacion de librerias

library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
library(gtable)


#Importar nuestra base de datos
list.files()
datos <- read.csv("all_genes_WT_vs_IrrKO.csv")
colnames(datos)

# Filtrar datos con padj no NA ni 0
datos <- datos %>%
  filter(!is.na(padj) & padj > 0)

#Construccion del grafico
# Agregar una columna a los datos para determinar el estado de la expresion
datos$Estado <- "Sin cambio"
# Si el log2Foldchange > 1 y qvalue < -log10(0.05) entonces "Sobreexpresado" 
datos$Estado[datos$log2FoldChange >= 1 & datos$padj < (0.05)] <- "Sobreexpresado"
# Si log2Foldchange < -1 y qvalue < -log10(0.05),} entonces "Subexpresado"
datos$Estado[datos$log2FoldChange <= -1 & datos$padj < (0.05)] <- "Subexpresado"

# Agregar el nombre de los genes
# Crear una columna "delabel" en la base de datos 
# que contenga el nombre de los genes expresados diferencialmente
datos$delabel <- NA
datos$delabel[datos$Estado != "NO"] <- datos$symbol[datos$Estado != "NO"]


# Insertar logarimos en las variables
# X/Y mejora la visualizacion
ggplot(data=datos, aes(x=datos$log2FoldChange, y=datos$padj)) + geom_point() #base

p = ggplot(data=datos, aes(x=log2FoldChange, y=-log10(padj), 
                           col=Estado)) + geom_point()+
  theme_classic(base_family = "Arial") + guides(color="none")
  #theme_classic()+ guides(color="none")
#labs(title = "Expresi?n diferencial entre GM-CSF vs 3CM178")+
#theme (plot.title = element_text (hjust = 0.5 ))+
#theme(legend.position="bottom")
p

# Agregar marco de referencia 
#de la expresion diferencial log2FoldChange y -log10qval 
#El log base 2 de 2 es 1...usar ese valor en lugar de o.6?
#El -log10(0.05)= 1.30
p2 = p + geom_vline(xintercept=c(-1, 1), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey")
p2

# Cambiar el color de los puntos
mycolors <- c("#1A85FF", "black", "#D41159")
names(mycolors) <- c("Subexpresado", "Sin cambio", "Sobreexpresado")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3






