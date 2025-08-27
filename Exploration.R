options(install.packages.check.source = "no")
BiocManager::install(ask = FALSE, update = FALSE)
# -- Libraries --

library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicRanges)
library(DESeq2)
library(tximport)
library(ggfortify)
library(AnnotationHub)
library(ensembldb)
library(enrichplot)
library(maser)
library(rtracklayer)
library(pheatmap)
library(Biostrings)
library(Rsamtools)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome)



# -- DATA --

# rMATS
se <- read_tsv("results/rmats/SE.MATS.JCEC.txt")
a5ss <- read_tsv("results/rmats/A5SS.MATS.JCEC.txt")
a3ss <- read_tsv("results/rmats/A3SS.MATS.JCEC.txt")
ri <- read_tsv("results/rmats/RI.MATS.JCEC.txt")
mxe  <- read_tsv("results/rmats/MXE.MATS.JCEC.txt")

# Salmon
metadata <- read.csv("samplesheet.csv")
files <- file.path("results/rnasplice/salmon", metadata$sample, "quant.sf")
names(files) <- metadata$sample
tx2gene <- read_tsv("tx2gene.tsv")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Filter by Significance --
sig_se <- se %>%
  dplyr::filter(FDR < 0.05 & abs(IncLevelDifference > 0.1))
sig_ri <- ri %>%
  dplyr::filter(FDR < 0.05 & abs(IncLevelDifference > 0.1))
sig_a3ss <- a3ss %>%
  dplyr::filter(FDR < 0.05 & abs(IncLevelDifference > 0.1))
sig_a5ss <- a5ss %>%
  dplyr::filter(FDR < 0.05 & abs(IncLevelDifference > 0.1))
sig_mxe <- mxe %>%
  dplyr::filter(FDR < 0.05 & abs(IncLevelDifference > 0.1))

table(c(
  rep("SE", nrow(sig_se)),
  rep("A5SS", nrow(sig_a5ss)),
  rep("A3SS", nrow(sig_a3ss)),
  rep("MXE", nrow(sig_mxe)),
  rep("RI", nrow(sig_ri))
))

# Parse IncLevels
split_psis <- function(x) {
  m <- strsplit(x, ",")
  maxlen <- max(lengths(m))
  m <- lapply(m, function(v) as.numeric(v[seq_len(maxlen)]))
  do.call(rbind, m)
}

# Set thresholds
FDR_CUTOFF  <- 0.05
DPSI_CUTOFF <- 0.1

# Clean group names from your metadata (assumes columns: sample, condition)
stopifnot(all(names(files) %in% metadata$sample))
group_lut <- setNames(metadata$condition, metadata$sample)
# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Visualization --

# PSI Distribution

all_sig <- bind_rows(
  sig_se %>% mutate(type="SE"),
  sig_a5ss %>% mutate(type="A5SS"),
  sig_a3ss %>% mutate(type="A3SS"),
  sig_ri %>% mutate(type="RI"),
  sig_mxe %>% mutate(type="MXE")
)

ggplot(all_sig, aes(x=IncLevelDifference, fill=type)) +
  geom_histogram(binwidth=0.05, alpha=0.7, position="identity") +
  theme_minimal() +
  labs(title="ΔPSI Distribution (significant events)",
       x="ΔPSI (Treatment - Control)",
       y="Count")

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Visualize rMATS --

counts_all <- c(
  SE   = nrow(se),
  A5SS = nrow(a5ss),
  A3SS = nrow(a3ss),
  MXE  = nrow(mxe),
  RI   = nrow(ri)
)

number_all <- c(SE =76719, RI = 7624, A5SS = 6880, A3SS = 10201, MXE = 13534)
label_all <- paste(names(number_all), number_all, sep="\n")

counts_sig <- c(
  SE   = nrow(sig_se),
  A5SS = nrow(sig_a5ss),
  A3SS = nrow(sig_a3ss),
  MXE  = nrow(sig_mxe),
  RI   = nrow(sig_ri)
)
number_sig <- c(SE =429, RI = 77, A5SS = 61, A3SS = 85, MXE = 58)
label_sig <- paste(names(number_sig), number_sig, sep="\n")
# Pie-chart of splicing methods observed
cols <- c(SE="#05668D", RI="#427AA1", A5SS="#EBF2FA", A3SS="#679436", MXE="#A5BE00")
op <- par(mfrow = c(1,2), mar = c(1,1,3,1))


pie(number_all, label_all,
    main = "All events", col = cols[names(number_all)], cex = 1.5)
pie(number_sig, label_sig,
    main = "Significant events", col = cols[names(counts_sig)], cex = 1.5)
par(op)


# Volcano-plot (Control vs. Treatment; SE)
volcano <- se %>%
  mutate(
    sig = case_when(
      FDR < 0.05 & IncLevelDifference > 0.1 ~ "Higher in Treatment",
      FDR < 0.05 & IncLevelDifference < -0.1 ~ "Higher in Control",
      TRUE ~ "Not significant"
    )
  )

ggplot(volcano, aes(x=IncLevelDifference, y=-log10(FDR), color=sig)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("red","blue","grey")) +
  geom_vline(xintercept=c(-0.1,0.1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal() +
  labs(title="ΔPSI Volcano Plot", x="ΔPSI (Treatment - Control)", y="-log10(FDR)")

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Significant instances vs. All --
bind_sig <- function(df, type) {
  df %>% mutate(type = type,
                sig = FDR < FDR_CUTOFF & abs(IncLevelDifference) > DPSI_CUTOFF)
}
all_events <- bind_rows(
  bind_sig(se,   "SE"),
  bind_sig(a5ss, "A5SS"),
  bind_sig(a3ss, "A3SS"),
  bind_sig(ri,   "RI"),
  bind_sig(mxe,  "MXE")
)
event_counts <- all_events %>%
  dplyr::count(type, name = "total") %>%
  left_join(all_events %>% dplyr::filter(sig) %>% dplyr::count(type, name = "significant"),
            by = "type") %>%
  mutate(significant = replace_na(significant, 0L)) %>%
  arrange(desc(total))
print(event_counts)
ggplot(event_counts, aes(x = reorder(type, -significant), y = significant, fill = type)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = significant), vjust = -0.25, size = 3.5) +
  labs(title = "Significant splicing events per type",
       y = "Count (FDR<0.05 & |ΔPSI|>0.05)", x = NULL) +
  theme_minimal(base_size = 12)

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- rMATS Heatmap --

# SE
se_sig <- se %>% dplyr::filter(FDR < FDR_CUTOFF, abs(IncLevelDifference) > DPSI_CUTOFF)

psi1 <- split_psis(se_sig$IncLevel1)
psi2 <- split_psis(se_sig$IncLevel2)

grp1_samples <- metadata$sample[metadata$condition == unique(metadata$condition)[1]]
grp2_samples <- metadata$sample[metadata$condition == unique(metadata$condition)[2]]
colnames(psi1) <- grp1_samples[seq_len(ncol(psi1))]
colnames(psi2) <- grp2_samples[seq_len(ncol(psi2))]
psi <- cbind(psi1, psi2)
rownames(psi) <- paste(se_sig$GeneID, se_sig$GeneID, sep = "|")

psi_z <- t(scale(t(psi), center = TRUE, scale = TRUE))

ann_col <- metadata %>%
  dplyr::select(sample, condition) %>%
  tibble::column_to_rownames("sample")
ann_colors <- list(
  condition = setNames(c("#06d6a0", "#118ab2"), unique(metadata$condition))
)

psi_z <- psi_z[, rownames(ann_col), drop = FALSE]
pheatmap(psi_z,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main = "Significant SE events: PSI (z-scored by event)")
# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Visualize TOP SE Events --
topN <- 6
picked <- se_sig %>% arrange(FDR, desc(abs(IncLevelDifference))) %>% slice_head(n = topN)
long_df <- map_dfr(seq_len(nrow(picked)), function(i) {
  r <- picked[i, ]
  s1 <- as.numeric(unlist(strsplit(r$IncLevel1, ",")))
  s2 <- as.numeric(unlist(strsplit(r$IncLevel2, ",")))
  tibble(
    Event = paste0(r$geneSymbol, " Exon ", r$exonID),
    Group = c(rep("Group1", length(s1)), rep("Group2", length(s2))),
    PSI   = c(s1, s2),
    FDR   = r$FDR
  )
})

lab_map <- c("Group1" = unique(metadata$condition)[1],
             "Group2" = unique(metadata$condition)[2])
long_df$Group <- lab_map[long_df$Group]
ggplot(long_df, aes(Group, PSI * 100)) +
  geom_boxplot(outlier.size = 0.8, width = 0.55) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.7) +
  facet_wrap(~ Event, scales = "free_y", ncol = 3) +
  labs(y = "PSI (%)", x = NULL,
       title = "Highlighted SE events (PSI per sample)") +
  theme_bw(base_size = 12)


# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- GO Enrichment --
# Dot-plot
sig_genes <- unique(c(sig_se$geneSymbol,
                      sig_a5ss$geneSymbol,
                      sig_a3ss$geneSymbol,
                      sig_mxe$geneSymbol,
                      sig_ri$geneSymbol))

ego <- enrichGO(sig_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05) %>%
  simplify()

enrichplot::dotplot(ego)
cnetplot(ego, showCategory = 2)

gene_description <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys=sig_genes,
                                          keytype="SYMBOL",
                                          columns=c("SYMBOL","GENENAME","ENTREZID"))
head(gene_description)


# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Differential Expression [DESeq2] --

sample_info <- metadata %>%
  mutate(condition = as.factor(condition))

## Create a Model
status_model <- as.formula(~ condition)

model.matrix(status_model, data = sample_info)

sample_info <- mutate(sample_info, Status = fct_relevel(condition, "HC"))
model.matrix(status_model, data = sample_info)

## Create a DESeq2 Object
dds_raw <- DESeqDataSetFromTximport(txi = txi,
                                    colData = sample_info,
                                    design = status_model)

## Filteration
keep <- rowSums(counts(dds_raw)) > 5
dds <- dds_raw[keep, ]

## DESeq Analysis
dds <- DESeq(dds)

# Results
res_status <- results(dds, alpha = 0.05)
res_status

res_condition_MDD_vs_HC <- res_status
rm(res_status)

# Identify Top Genes
topGenesIvU <- res_condition_MDD_vs_HC %>%
  as.data.frame() %>%
  rownames_to_column("GeneID") %>%
  top_n(100, wt = -padj)

head(topGenesIvU)

sum(res_condition_MDD_vs_HC$padj < 0.05, na.rm = TRUE)

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Database Query for DESeq2 --

ah <- AnnotationHub()

ensdb <- query(ah, c("EnsDb", "Homo sapien", "102"))[[1]]

annot <- genes(ensdb, return.type = "data.frame")
colnames(annot)

annot <- annot %>%
  dplyr::select(gene_id, gene_name, entrezid) %>%
  dplyr::filter(gene_id %in% rownames(res_condition_MDD_vs_HC))

# Apply to our Data


# ----------------------------------------------------------------
# ----------------------------------------------------------------

# -- Motif Matching --

