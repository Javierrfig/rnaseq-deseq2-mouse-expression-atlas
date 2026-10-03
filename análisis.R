# ==============================================================================
# RNA-Seq differential expression with DESeq2 - GXA E-MTAB-10413 (Mus musculus)
# Knockout SLC38A10 vs wild type, en 4 condiciones (estres ambiental + tiempo)
# Adaptado de un tutorial de RNA-seq con DESeq2 (OMGenomics)
# ==============================================================================

# ---- Setup -------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2", "EnhancedVolcano", "biomaRt")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(biomaRt)
library(patchwork)

out_dir <- "~/rna-seq"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Umbrales y genes de interes (cambia aqui y se aplican en todo el script) -
padj_cut <- 0.05
lfc_cut  <- 7
genes_check <- c("Arsi", "Col8a1", "Madcam1", "Anxa1", "A730049H05Rik")
gxa_genes   <- c(genes_check, "Glycam1", "Acta2")

# ---- Download data -----------------------------------------------------------
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10413/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
meta_raw = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10413/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(counts)
head(meta_raw)

# ---- Wrangle counts ----------------------------------------------------------
rownames(counts) = counts$Gene.ID
genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]

# ---- Wrangle metadata --------------------------------------------------------
meta_raw = meta_raw[!duplicated(meta_raw$Run), ]
rownames(meta_raw) = meta_raw$Run

find_col = function(df, pattern) {
  hits = grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)
  hits[!grepl("ontology|unit", hits, ignore.case = TRUE)]
}
geno_col   = find_col(meta_raw, "genotype")
stress_col = find_col(meta_raw, "stress")
time_col   = find_col(meta_raw, "time")
print(list(genotype = geno_col, stress = stress_col, time = time_col))
if (length(geno_col) < 1 || length(stress_col) < 1 || length(time_col) < 1) {
  print(colnames(meta_raw))
  stop("No encontre las columnas de genotype/estres/tiempo. Revisa los nombres impresos arriba.")
}

g = trimws(as.character(meta_raw[[geno_col[1]]]))
g = ifelse(grepl("wild", g, ignore.case = TRUE), "wildtype", make.names(g))
lv = c(intersect("wildtype", g), sort(setdiff(unique(g), "wildtype")))

metadata = data.frame(
  genotype  = factor(g, levels = lv),
  condition = paste(trimws(meta_raw[[stress_col[1]]]), "at", trimws(meta_raw[[time_col[1]]])),
  row.names = rownames(meta_raw)
)
print(table(metadata$genotype, metadata$condition, useNA = "ifany"))
stopifnot("Quedan NA en genotype" = !anyNA(metadata$genotype))

common = intersect(colnames(counts), rownames(metadata))
stopifnot("Ninguna muestra coincide entre counts y metadata" = length(common) > 0)
counts = counts[, common]
metadata = metadata[common, , drop = FALSE]

ref  = levels(metadata$genotype)[1]
case = levels(metadata$genotype)[2]
message("Contraste: ", case, " vs ", ref)

# ---- Spot check del gen knockout ---------------------------------------------
gene_name = "Slc38a10"
gene_id = genes$Gene.ID[genes$Gene.Name == gene_name]
if (length(gene_id) == 1) {
  gene_data = cbind(metadata, counts = as.numeric(counts[gene_id, ]))
  print(ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.1, size = 2) +
          facet_wrap(~condition) + ggtitle(gene_name))
} else message("Gen ", gene_name, " no encontrado; se omite el spot check.")

# ---- Analisis A: todas las condiciones juntas (~genotype) --------------------
# Vista general; NO es comparable con GXA (mezcla las 4 condiciones)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~genotype)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res = results(dds, contrast = c("genotype", case, ref), alpha = padj_cut)

res_df = as.data.frame(res)
res_df$Gene.ID = rownames(res_df)
res_df = merge(res_df, genes, by = "Gene.ID")
top_genes_A = head(res_df$Gene.Name[order(res_df$padj)], 4)
print(res_df[res_df$Gene.Name %in% top_genes_A, ])

# ---- Analisis B: una comparacion por condicion (como GXA) --------------------
res_by_cond = list()
for (cn in unique(metadata$condition)) {
  idx = metadata$condition == cn
  md  = droplevels(metadata[idx, , drop = FALSE])
  if (nlevels(md$genotype) < 2 || any(table(md$genotype) < 2)) {
    message("Se omite '", cn, "': faltan muestras de algun genotype.")
    next
  }
  d <- DESeqDataSetFromMatrix(countData = counts[, idx], colData = md, design = ~genotype)
  d <- d[rowSums(counts(d)) > 10, ]
  d <- DESeq(d)
  res_by_cond[[cn]] = results(d, contrast = c("genotype", case, ref), alpha = padj_cut)
}

# Tabla larga con todas las condiciones (la tabla principal)
res_cond_df = do.call(rbind, lapply(names(res_by_cond), function(cn) {
  r = as.data.frame(res_by_cond[[cn]])
  r$Gene.ID = rownames(r)
  r = merge(r, genes, by = "Gene.ID")
  r$condition = cn
  r
}))

# Genes de la captura de GXA, por condicion
print(res_cond_df[res_cond_df$Gene.Name %in% gxa_genes,
                  c("condition", "Gene.Name", "log2FoldChange", "padj")])

# Filtro tipo GXA (|log2FC| >= lfc_cut y p ajustado <= padj_cut)
gxa_like = res_cond_df[!is.na(res_cond_df$padj) &
                         res_cond_df$padj <= padj_cut &
                         abs(res_cond_df$log2FoldChange) >= lfc_cut, ]
gxa_like = gxa_like[order(gxa_like$condition, gxa_like$log2FoldChange), ]
print(gxa_like[, c("condition", "Gene.Name", "log2FoldChange", "padj")])

# Tabla ancha tipo GXA: un gen por fila, una condicion por columna
tabla_gxa = reshape(gxa_like[, c("Gene.ID", "Gene.Name", "condition", "log2FoldChange")],
                    idvar = c("Gene.ID", "Gene.Name"), timevar = "condition",
                    direction = "wide")
names(tabla_gxa) = sub("log2FoldChange\\.", "", names(tabla_gxa))
tabla_gxa$n_cond = rowSums(!is.na(tabla_gxa[, -(1:2)]))
tabla_gxa = tabla_gxa[order(tabla_gxa$n_cond), ]
rownames(tabla_gxa) = NULL

if (interactive()) {
  View(res_cond_df[res_cond_df$Gene.Name %in% gxa_genes,
                   c("condition", "Gene.Name", "log2FoldChange", "padj")])
  View(gxa_like[, c("condition", "Gene.Name", "log2FoldChange", "padj")])
  View(tabla_gxa)
}

# ---- Visualizacion -----------------------------------------------------------
# Figura 1: MA plot + volcano (analisis A)
ma_data <- plotMA(res, returnData = TRUE)
p_ma <- ggplot(ma_data, aes(mean, lfc, color = isDE)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c("grey60", "red")) +
  labs(title = "MA plot", x = "Media normalizada", y = "log2 fold change") +
  theme_minimal() + theme(legend.position = "none")
p_volcano <- EnhancedVolcano(res, lab = rownames(res),
                             x = "log2FoldChange", y = "pvalue")
combinado <- p_ma + p_volcano
print(combinado)
ggsave(file.path(fig_dir, "ma_volcano.png"), combinado, width = 14, height = 6, dpi = 300)

# Figura 2: un volcano por condicion (analisis B)
volcanos = lapply(names(res_by_cond), function(cn)
  EnhancedVolcano(res_by_cond[[cn]], lab = rep("", nrow(res_by_cond[[cn]])),
                  x = "log2FoldChange", y = "pvalue", title = cn,
                  subtitle = NULL, caption = NULL))
por_condicion <- wrap_plots(volcanos, ncol = 2)
print(por_condicion)
ggsave(file.path(fig_dir, "volcano_por_condicion.png"), por_condicion,
       width = 14, height = 12, dpi = 300)

# ---- Guardar resultados -----------------------------------------------------
write.csv(res_cond_df, file.path(out_dir, "deseq_por_condicion.csv"), row.names = FALSE)
write.csv(tabla_gxa, file.path(out_dir, "tabla_tipo_gxa.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# ---- Coordenadas de los genes (Ensembl, Mus musculus) ------------------------
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

# Alternativa local si Ensembl no responde (descomenta y ejecuta; omite la linea
# colnames(all.genes)[1] de abajo):
# BiocManager::install("EnsDb.Mmusculus.v79")
# library(EnsDb.Mmusculus.v79)
# g2 <- as.data.frame(genes(EnsDb.Mmusculus.v79))
# all.genes <- data.frame(Gene.ID = g2$gene_id, chromosome_name = g2$seq_name,
#                         start_position = g2$start, end_position = g2$end)

stopifnot("Ensembl no responde; usa la alternativa local" = !is.null(all.genes))
colnames(all.genes)[1] = "Gene.ID"
message("Genes con coordenadas: ", round(100 * mean(res_cond_df$Gene.ID %in% all.genes$Gene.ID), 1), "%")

all.genes <- all.genes[all.genes$chromosome_name %in% c(1:19, "X", "Y"), ]
all.genes$chromosome_name <- paste0("chr", all.genes$chromosome_name)

# Resultados con coordenadas
merged_data <- merge(all.genes, res_cond_df, by = "Gene.ID")
merged_data <- merged_data[!is.na(merged_data$log2FoldChange), ]
message("Filas en merged_data: ", nrow(merged_data))

# Analisis A con coordenadas
merged_A <- merge(all.genes, res_df, by = "Gene.ID")
write.csv(merged_A, file.path(out_dir, "deseq.csv"), row.names = FALSE)
write.csv(merged_A[merged_A$Gene.Name %in% top_genes_A, ],
          file.path(out_dir, "deseq_subset.csv"), row.names = FALSE)

# ---- Genes vecinos del gen knockout: hay mas alterados de lo esperable? ------
# H0: los genes cercanos a Slc38a10 no se alteran mas que el resto del genoma
target    <- "Slc38a10"
window_bp <- 1e6
lfc_min   <- 1

tg <- merged_data[merged_data$Gene.Name == target, ][1, ]
stopifnot("No encontre el gen en merged_data" = !is.na(tg$Gene.ID))
chr_t <- tg$chromosome_name
pos_t <- (tg$start_position + tg$end_position) / 2
message(target, " esta en ", chr_t, ", posicion ~", round(pos_t / 1e6, 2), " Mb")

neighbor_test <- function(lfc) {
  out <- lapply(unique(merged_data$condition), function(cn) {
    d <- merged_data[merged_data$condition == cn & !is.na(merged_data$padj) &
                       merged_data$Gene.Name != target, ]
    d$pos  <- (d$start_position + d$end_position) / 2
    d$near <- d$chromosome_name == chr_t & abs(d$pos - pos_t) <= window_bp
    d$sig  <- d$padj <= padj_cut & abs(d$log2FoldChange) >= lfc
    if (sum(d$near) < 5) { message("Pocos genes en la ventana: ", cn); return(NULL) }
    tab <- table(near = factor(d$near, levels = c(TRUE, FALSE)),
                 sig  = factor(d$sig,  levels = c(TRUE, FALSE)))
    ft <- fisher.test(tab)
    data.frame(lfc_min = lfc, condition = cn,
               genes_cerca = sum(d$near),
               alterados_cerca = sum(d$near & d$sig),
               pct_cerca = round(100 * mean(d$sig[d$near]), 1),
               pct_resto = round(100 * mean(d$sig[!d$near]), 1),
               odds_ratio = unname(ft$estimate),
               p_fisher = ft$p.value)
  })
  res <- do.call(rbind, out)
  if (!is.null(res)) res$p_ajustado <- p.adjust(res$p_fisher, method = "BH")
  res
}

local_res <- neighbor_test(lfc_min)
if (!is.null(local_res)) {
  print(local_res)
  write.csv(local_res, file.path(out_dir, "vecinos_gen_knockout.csv"), row.names = FALSE)
}

# ---- Analisis de sensibilidad ------------------------------------------------
lfc_grid <- c(0, 0.68, 1, 2)
sens_res <- do.call(rbind, lapply(lfc_grid, neighbor_test))
if (!is.null(sens_res)) {
  print(sens_res[, c("lfc_min", "condition", "genes_cerca", "alterados_cerca",
                     "pct_cerca", "pct_resto", "odds_ratio", "p_ajustado")])
  write.csv(sens_res, file.path(out_dir, "vecinos_sensibilidad.csv"), row.names = FALSE)

  print(aggregate(p_ajustado ~ lfc_min, data = sens_res,
                  FUN = function(x) sum(x < 0.05, na.rm = TRUE)))

  sens_plot <- sens_res[is.finite(sens_res$odds_ratio) & sens_res$odds_ratio > 0, ]
  sens_plot$significativo <- sens_plot$p_ajustado < 0.05
  p_sens <- ggplot(sens_plot, aes(factor(lfc_min), odds_ratio, color = condition,
                                  shape = significativo, group = condition)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line() + geom_point(size = 3) +
    scale_y_log10() +
    labs(title = paste("Sensibilidad: vecinos de", target),
         x = "|log2FC| minimo para llamar 'alterado'",
         y = "Odds ratio (escala log)", shape = "p ajustado < 0.05")
  print(p_sens)
  ggsave(file.path(fig_dir, "sensibilidad_vecinos.png"), p_sens, width = 9, height = 6, dpi = 300)
}

# Grafico del cromosoma
plot_df <- merged_data[merged_data$chromosome_name == chr_t & !is.na(merged_data$padj), ]
plot_df$pos_mb <- (plot_df$start_position + plot_df$end_position) / 2 / 1e6
plot_df$sig <- plot_df$padj <= padj_cut & abs(plot_df$log2FoldChange) >= lfc_min
p_local <- ggplot(plot_df, aes(pos_mb, log2FoldChange, color = sig)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = pos_t / 1e6, linetype = "dashed") +
  scale_color_manual(values = c("grey60", "red")) +
  facet_wrap(~condition) +
  labs(title = paste("Genes en", chr_t), x = "Posicion (Mb)", y = "log2FC", color = "alterado")
print(p_local)
ggsave(file.path(fig_dir, "cromosoma_gen_knockout.png"), p_local, width = 12, height = 7, dpi = 300)
