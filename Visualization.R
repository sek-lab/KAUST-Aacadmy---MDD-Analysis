# ---- packages ----
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(ggplot2); library(tidyr); library(purrr)
})

# ---- load rMATS (your existing files) ----
se   <- read_tsv("results/rmats/SE.MATS.JCEC.txt",  show_col_types = FALSE)  %>% mutate(event="SE")
a5ss <- read_tsv("results/rmats/A5SS.MATS.JCEC.txt",show_col_types = FALSE)  %>% mutate(event="A5SS")
a3ss <- read_tsv("results/rmats/A3SS.MATS.JCEC.txt",show_col_types = FALSE)  %>% mutate(event="A3SS")
ri   <- read_tsv("results/rmats/RI.MATS.JCEC.txt",  show_col_types = FALSE)  %>% mutate(event="RI")
mxe  <- read_tsv("results/rmats/MXE.MATS.JCEC.txt", show_col_types = FALSE)  %>% mutate(event="MXE")

rmats <- bind_rows(se,a5ss,a3ss,ri,mxe)

# ---- quick counts by event type ----
counts_all <- rmats %>% dplyr::count(event, name="N_all")
counts_sig <- rmats %>% dplyr::filter(FDR <= 0.05 & (abs(IncLevelDifference) > 0.1)) %>% dplyr::count(event, name="N_sig")
left_join(counts_all, counts_sig, by="event") %>% replace_na(list(N_sig=0))
# (prints a small table)
rmats <- rmats %>%
  dplyr::mutate(nlogFDR = -log10(pmax(FDR, 1e-300))) %>%
  dplyr::filter(!(event == "RI" & nlogFDR > 50))   # remove RI points with -log10(FDR) > 50

# ---- volcano (ΔPSI vs -log10FDR), faceted ----
rmats %>%
  mutate(
    nlogFDR = -log10(pmax(FDR, 1e-300)),
    sig = FDR < 0.05 & abs(IncLevelDifference) >= 0.1,   # mark significant
    direction = case_when(
      sig & IncLevelDifference > 0 ~ "MDD higher inclusion",
      sig & IncLevelDifference < 0 ~ "Control higher inclusion",
      TRUE ~ "Not Significant"
    )
  ) %>%
  ggplot(aes(x = IncLevelDifference, y = nlogFDR, color = direction)) +
  geom_point(alpha=.6, size=1) +
  geom_vline(xintercept = c(-0.2,0.2), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  facet_wrap(~event, nrow=2, scales="free_y") +
  scale_color_manual(values=c(
    "MDD higher inclusion" = "#A5BE00",
    "Control higher inclusion" = "#377eb8",
    "Not Significant" = "grey70"
  )) +
  labs(
    title = "ΔPSI volcano by event type",
    x = "ΔPSI (MDD – Control)",
    y = "-log10(FDR)",
    color = "Direction"
  ) +
  theme_minimal(base_size = 12)


ggsave("volcano_rmats_facets.pdf", width=10, height=6)

# ---- helper: parse per-sample PSI strings to numeric vector ----
suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(stringr); library(readr)
})

# helpers
split_num <- function(x) {
  if (is.na(x) || x == "") return(numeric())
  suppressWarnings(as.numeric(str_split(x, "[,;]")[[1]]))
}
sum_counts <- function(x) {
  v <- split_num(x); v[is.na(v)] <- 0; sum(v)
}

candidates <- rmats %>%
  # vectorized parts
  dplyr::mutate(
    meanPSI_1 = map_dbl(IncLevel1, ~mean(split_num(.x), na.rm = TRUE)),
    meanPSI_2 = map_dbl(IncLevel2, ~mean(split_num(.x), na.rm = TRUE))
  ) %>%
  # row-wise for counts & consistency
  rowwise() %>%
  dplyr::mutate(
    reads_g1    = sum_counts(IJC_SAMPLE_1) + sum_counts(SJC_SAMPLE_1),
    reads_g2    = sum_counts(IJC_SAMPLE_2) + sum_counts(SJC_SAMPLE_2),
    reads_total = reads_g1 + reads_g2,
    frac_high_1 = mean(split_num(IncLevel1) > 0.5, na.rm = TRUE),
    frac_high_2 = mean(split_num(IncLevel2) > 0.5, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # biomarker-like filter (tune thresholds if you want)
  dplyr::filter(
    FDR < 0.05,
    abs(IncLevelDifference) >= 0.10,
    reads_total >= 50,
    abs(frac_high_2 - frac_high_1) >= 0.5
  ) %>%
  dplyr::arrange(desc(abs(IncLevelDifference)))

# outputs (same as your originals)
write_tsv(candidates, "rmats_candidates_all.tsv")

top20 <- candidates %>%
  group_by(event) %>%
  slice_max(order_by = abs(IncLevelDifference), n = 20, with_ties = FALSE) %>%
  ungroup()
write_tsv(top20, "rmats_top20_by_event.tsv")

gene_hits <- candidates %>%
  dplyr::count(geneSymbol, sort = TRUE, name = "n_events")
write_tsv(gene_hits, "rmats_candidate_gene_counts.tsv")






# ---- foreground / background genes ----
sig_tbl <- rmats %>% dplyr::filter(FDR < 0.05, abs(IncLevelDifference) >= 0.1)
sig_genes <- unique(na.omit(sig_tbl$geneSymbol))
bg_genes  <- unique(na.omit(rmats$geneSymbol))

# ---- hypergeometric against local GO text files (paper-style) ----
scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication",
                 manual_pal(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
substrRight <- function(x, n) substr(x, nchar(x)-n+1, nchar(x))

run_hypergeometric <- function(observed_gene_list, bg_gene_list, dir_database="."){
  db_vec <- c("GO_Biological_Process_2025","GO_Cellular_Component_2025","GO_Molecular_Function_2025")
  db_label <- c("Biological_Process","Cellular_Component","Molecular_Function"); names(db_label) <- db_vec
  out <- list()
  for (db in db_vec){
    con <- file(file.path(dir_database, paste0(db,".txt")), "r")
    while(length(oneLine <- readLines(con, n=1, warn=FALSE))>0){
      x <- strsplit(oneLine, "\t")[[1]]
      term  <- x[1]; genes <- x[-1]; genes <- genes[genes!=""]
      total <- toupper(bg_gene_list)
      g1 <- toupper(observed_gene_list)
      g2 <- intersect(toupper(genes), total)
      ov <- intersect(g1, g2)
      if (!length(ov)) next
      pval <- phyper(length(ov)-1, length(g2), length(total)-length(g2), length(g1), lower.tail=FALSE)
      odds <- (length(ov)/length(g1)) / (length(g2)/length(total))
      out[[length(out)+1]] <- data.frame(
        Term=term, P.value=pval, Odds.Ratio=odds,
        Genes=paste(ov, collapse=";"),
        Overlap_number=length(ov), Term_gene_number=length(g2),
        db=db_label[db], stringsAsFactors=FALSE)
    }
    close(con)
  }
  if (!length(out)) return(data.frame())
  df <- bind_rows(out) %>% mutate(Adjusted.P.value=p.adjust(P.value, "fdr"))
  df
}

df_enrich <- run_hypergeometric(sig_genes, bg_genes, ".")
stopifnot(nrow(df_enrich)>0)

# pick top per category (like the paper)
pick_top <- function(d, cutoff=0.05){
  d <- arrange(d, Adjusted.P.value)
  if (nrow(d)==0) return(d)
  if (min(d$Adjusted.P.value) <= cutoff){
    keep <- min(10, max(min(5, nrow(d)), sum(d$Adjusted.P.value <= cutoff)))
  } else keep <- min(5, nrow(d))
  d[seq_len(keep),]
}
plot_df <- bind_rows(
  df_enrich %>% dplyr::filter(db=="Biological_Process")  %>% pick_top(),
  df_enrich %>% dplyr::filter(db=="Cellular_Component")  %>% pick_top(),
  df_enrich %>% dplyr::filter(db=="Molecular_Function")  %>% pick_top()
)

# shorten long names
long <- nchar(plot_df$Term) > 100
plot_df$Term[long] <- paste(substr(plot_df$Term[long],1,100), substrRight(plot_df$Term[long],12), sep="… ")

p <- ggplot(plot_df,
            aes(x=reorder(Term, -log10(Adjusted.P.value)),
                y=-log10(Adjusted.P.value), fill=db, alpha=Odds.Ratio)) +
  geom_col(width=.6, position=position_dodge2(width=.7, preserve="single")) +
  facet_grid(rows=vars(db), scales="free_y", space="free_y", switch="both") +
  scale_fill_manual(
    values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c",
               "#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"),
    guide = "none"
  ) +
  scale_alpha_continuous(range=c(0.3,1.0)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_text(aes(label=paste0(Overlap_number,"/",Term_gene_number),
                y=-log10(Adjusted.P.value) + 0.05),
            size=3, color="grey30", hjust=0) +
  coord_flip() +
  labs(title="GO enrichment", x="Terms", y="-log10(adjusted p value)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(),
        strip.background=element_blank(),
        strip.text.y.left=element_text(angle=0),
        legend.position="bottom")

ggsave("plot_enrich_GO_faceted.pdf", p, width=10, height=5)
write_tsv(df_enrich %>% arrange(db, Adjusted.P.value), "go_enrich_full.tsv")

