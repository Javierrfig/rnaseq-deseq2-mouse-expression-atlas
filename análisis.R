# ==============================================================================
# RNA-Seq differential expression with DESeq2 - GXA E-MTAB-10413
# Version adaptada: detecta automaticamente los valores de genotype
# ==============================================================================

# ---- Setup -------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2", "EnhancedVolcano", "biomaRt")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(biomaRt)

dir.create("~/Documents/rna-seq", recursive = TRUE, showWarnings = FALSE)

# ---- Download data -----------------------------------------------------------
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10413/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10413/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(counts)
head(metadata)

# ---- Wrangle counts ----------------------------------------------------------
rownames(counts) = counts$Gene.ID
genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]

# ---- Wrangle metadata --------------------------------------------------------
# Una fila por muestra (Run)
metadata = metadata[!duplicated(metadata$Run), ]
rownames(metadata) = metadata$Run

# Encuentra la columna de genotype sin depender del nombre exacto
geno_col = grep("genotype", colnames(metadata), ignore.case = TRUE, value = TRUE)
print(geno_col)
stopifnot("No hay columna de genotype en metadata" = length(geno_col) >= 1)
metadata = metadata[, geno_col[1], drop = FALSE]
colnames(metadata) = "genotype"

# Mira los valores reales (revisalos en la consola)
print(table(metadata$genotype, useNA = "ifany"))

# Recodifica: cualquier "wild type" -> wildtype; el resto, nombre sin espacios
g = trimws(as.character(metadata$genotype))
g = ifelse(grepl("wild", g, ignore.case = TRUE), "wildtype", make.names(g))
if (!"wildtype" %in% g) message("Aviso: no hay 'wild type'; el nivel de referencia sera: ", sort(unique(g))[1])
lv = c(intersect("wildtype", g), sort(setdiff(unique(g), "wildtype")))
metadata$genotype = factor(g, levels = lv)
print(metadata)

stopifnot("Quedan NA en genotype" = !anyNA(metadata$genotype))

# Muestras de metadata y de counts deben coincidir y estar en el mismo orden
common = intersect(colnames(counts), rownames(metadata))
stopifnot("Ninguna muestra coincide entre counts y metadata" = length(common) > 0)
counts = counts[, common]
metadata = metadata[common, , drop = FALSE]

# ---- Spot check de un gen (cambia gene_name por tu gen de interes) -----------
gene_name = "SNAI1"
gene_id = genes$Gene.ID[genes$Gene.Name == gene_name]
if (length(gene_id) == 1) {
  gene_data = cbind(metadata, counts = as.numeric(counts[gene_id, ]))
  print(ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot())
} else message("Gen ", gene_name, " no encontrado; se omite el spot check.")

# ---- Run DESeq ---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~genotype)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Compara cada nivel contra wildtype (o el primer nivel)
ref = levels(metadata$genotype)[1]
case = levels(metadata$genotype)[2]
message("Contraste: ", case, " vs ", ref)
res = results(dds, contrast = c("genotype", case, ref), alpha = 1e-5)
res

# ---- Merge con nombres de genes ----------------------------------------------
res_df = as.data.frame(res)
res_df$Gene.ID = rownames(res_df)
res_df = merge(res_df, genes, by = "Gene.ID")
head(res_df)

# Genes a revisar: los 4 mas significativos
genes_to_check = head(res_df$Gene.Name[order(res_df$padj)], 4)
res_df[res_df$Gene.Name %in% genes_to_check, ]

# ---- Visualizacion -----------------------------------------------------------
plotMA(res)
EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "pvalue")

# ---- Coordenadas (Ensembl) ---------------------------------------------------
all.genes <- NULL
for (m in c("useast", "asia", "www")) {
  message("Probando mirror: ", m)
  all.genes <- tryCatch({
    ensembl <- useEnsembl(biomart = "genes",
                          dataset = "mmusculus_gene_ensembl",
                          mirror = m)
    getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                         "start_position", "end_position"),
          mart = ensembl)
  }, error = function(e) NULL)
  if (!is.null(all.genes)) break
}

# Alternativa local si Ensembl no responde (descomenta y ejecuta):
# BiocManager::install("EnsDb.Mmusculus.v79")
# library(EnsDb.Mmusculus.v79)
# g <- as.data.frame(genes(EnsDb.Mmusculus.v79))
# all.genes <- data.frame(Gene.ID = g$gene_id, chromosome_name = g$seq_name,
#                         start_position = g$start, end_position = g$end)
# (en ese caso omite la linea colnames(all.genes)[1] de abajo)

stopifnot("Ensembl no responde; usa la alternativa local" = !is.null(all.genes))
colnames(all.genes)[1] = "Gene.ID"

merged_data <- merge(all.genes, res_df, by = "Gene.ID")
merged_data$chromosome_name <- paste0("chr", merged_data$chromosome_name)
write.csv(merged_data, "~/Documents/rna-seq/deseq.csv", row.names = FALSE)

merged_data_subset = merged_data[merged_data$Gene.Name %in% genes_to_check, ]
write.csv(merged_data_subset, "~/Documents/rna-seq/deseq_subset.csv", row.names = FALSE)
