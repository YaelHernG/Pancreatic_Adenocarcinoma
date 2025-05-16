####################################
# Author: Paola Albarran Godoy, Yael Daniel Hernandez Gonzalez, Ariadna Badia Zamudio
# Date: 10 de mayo 2025
# Description:

#Input : 

#Genes drivers identification

rm(list = ls())
# Usa data.table para una carga rápida
library(data.table)
maf <- fread("/home/yael/2025-2/data_mutations.txt", data.table = FALSE)
# Extraer las columnas necesarias y renombrarlas
muts <- maf[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                "Reference_Allele", "Tumor_Seq_Allele2")]
colnames(muts) <- c("sampleID", "chr", "pos", "ref", "mut")
# Guardar como archivo de entrada para dndscv
write.table(muts, file = "pancreas.5col", sep = "\t", quote = FALSE, row.names = FALSE)

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("im3sanger/dndscv")

# Cargar el paquete
library(dndscv)

# Leer el archivo con las mutaciones
muts = read.table("pancreas.5col", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Número de muestras únicas
num_muestras <- length(unique(muts$sampleID))
cat("Número de muestras:", num_muestras, "\n")

# Número total de mutaciones
num_mutaciones <- nrow(muts)
cat("Número total de mutaciones:", num_mutaciones, "\n")

# Identificar hipermutadores (más de 500 mutaciones)
hipermutadores <- names(mut_por_muestra[mut_por_muestra > 500])
cat("Muestras hipermutadoras:", hipermutadores, "\n")

# Ejecutar dndscv con parámetros estándar
dout <- dndscv(muts, 
               max_muts_per_gene_per_sample = 3,
               max_coding_muts_per_sample = 500,
               outmats = TRUE)

# Ver los nombres de los objetos en el resultado
names(dout)

# Mostrar genes con evidencia de selección positiva (q < 0.1)
genes_significativos <- dout$sel_cv[which(dout$sel_cv$qglobal_cv < 0.1), ]
print(genes_significativos)

# Ver las estimaciones globales de dN/dS
print(dout$globaldnds)


#Genes under positive selection
#TP53
# Coeficiente de selección para mutaciones missense en TP53
wmis_TP53 = dout$sel_cv$wmis_cv[dout$sel_cv$gene_name == "TP53"]
n_mis_TP53 = dout$sel_cv$n_mis[dout$sel_cv$gene_name == "TP53"]

# Proporción de mutaciones seleccionadas
selected_prop = (wmis_TP53 - 1) / wmis_TP53 
selected_muts = n_mis_TP53 * selected_prop
cat("Mutaciones missense bajo selección positiva en TP53:", round(selected_muts), "\n")

# Confianza del coeficiente
geneci(dout, gene_list = "TP53")

# Mutaciones observadas vs esperadas en SMAD4
dout$genemuts[which(dout$genemuts$gene_name == "TP53"), ]

#SMAD4
# Coeficiente de selección para mutaciones missense en SMAD4
wmis_SMAD4 = dout$sel_cv$wmis_cv[dout$sel_cv$gene_name == "SMAD4"]
n_mis_SMAD4 = dout$sel_cv$n_mis[dout$sel_cv$gene_name == "SMAD4"]

# Proporción de mutaciones seleccionadas
selected_prop = (wmis_SMAD4 - 1) / wmis_SMAD4 # El 76% de las mutaciones son "exceso", probablemente conductoras
selected_muts = n_mis_SMAD4 * selected_prop
cat("Mutaciones missense bajo selección positiva en SMAD4:", round(selected_muts), "\n")

# Confianza del coeficiente
geneci(dout, gene_list = "SMAD4")

# Mutaciones observadas vs esperadas en SMAD4
dout$genemuts[which(dout$genemuts$gene_name == "SMAD4"), ]

#KRAS
# Coeficiente de selección para mutaciones missense en SMAD4
wmis_KRAS = dout$sel_cv$wmis_cv[dout$sel_cv$gene_name == "KRAS"]
n_mis_KRAS = dout$sel_cv$n_mis[dout$sel_cv$gene_name == "KRAS"]

# Proporción de mutaciones seleccionadas
selected_prop = (wmis_KRAS - 1) / wmis_KRAS
selected_muts = n_mis_KRAS * selected_prop
cat("Mutaciones missense bajo selección positiva en KRAS:", round(selected_muts), "\n")

# Confianza del coeficiente
geneci(dout, gene_list = "KRAS")

# Mutaciones observadas vs esperadas en SMAD4
dout$genemuts[which(dout$genemuts$gene_name == "KRAS"), ]

##########
#Visualization
# Cargar el paquete maftools
library(maftools)
# Leer el archivo MAF
pancreas_maf <- read.maf(maf = "/home/yael/2025-2/data_mutations.txt")

oncoplot(pancreas_maf, genes = c("TP53","KRAS","SMAD4","CDKN2A","MUC16"), draw_titv
         = TRUE, showTumorSampleBarcodes=FALSE)

lollipopPlot(maf=pancreas_maf, gene="TP53", showMutationRate=TRUE)
lollipopPlot(maf=pancreas_maf, gene="KRAS", showMutationRate=TRUE)
lollipopPlot(maf=pancreas_maf, gene="SMAD4", showMutationRate=TRUE)
lollipopPlot(maf=pancreas_maf, gene="CDKN2A", showMutationRate=TRUE)
lollipopPlot(maf=pancreas_maf, gene="MUC16", showMutationRate=TRUE)


#Signalign pathways
# Instalar y cargar pathview
BiocManager::install("pathview")
BiocManager::install("org.Hs.eg.db")
library(pathview)
library(org.Hs.eg.db)
# Definir los valores para cada gen
gene.values <- c("KRAS" = 66,    # Porcentaje de muestras alteradas
                 "TP53" = 61,
                 "SMAD4" = 21,
                 "CDKN2A" = 20,
                 "MUC16" = 7)

# Mapeo de los nombres a ENTREZ IDs usando org.Hs.eg.db
entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = names(gene.values),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Asignar los valores a los IDs
names(gene.values) <- entrez_ids

# Generar la visualización para las vías de señalización
pathview(gene.data = gene.values, 
         pathway.id = "hsa05212", 
         species = "hsa", 
         gene.idtype = "entrez",
         kegg.native = TRUE, 
         limit = list(gene = c(0, 70), cpd = 1),
         low = "lightgreen", mid = "yellow", high = "red")