# === deps ===
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(stringr)
  library(tidyr); library(glmnet); library(pROC)
})

# ======================
# USER INPUTS (edit)
# ======================
paths <- list(
  SE   = "results/rmats/SE.MATS.JCEC.txt",
  MXE  = "results/rmats/MXE.MATS.JCEC.txt",
  A5SS = "results/rmats/A5SS.MATS.JCEC.txt",
  A3SS = "results/rmats/A3SS.MATS.JCEC.txt",
  RI   = "results/rmats/RI.MATS.JCEC.txt"
)

# Filters for “candidate” events
FDR_MAX   <- 0.10    # relax/tighten as needed
DPSI_MIN  <- 0.15
READS_MIN <- 10

# Modeling knobs
TOP_UNI       <- 10000   # keep top-K univariate features before LASSO
MAX_NA_PROP   <- 0.50  # drop features with >50% missing PSI across samples
USE_ALPHA     <- 1.0   # start with LASSO; we’ll fallback to 0.5 if empty
N_FOLDS       <- 3

# Outputs
OUT_DIR   <- "biomarker_combo"
dir.create(OUT_DIR, showWarnings = FALSE)
FN_PANEL  <- file.path(OUT_DIR, "selected_panel.tsv")
FN_PSI    <- file.path(OUT_DIR, "psi_matrix.tsv")
FN_META   <- file.path(OUT_DIR, "metadata.tsv")
FN_UNI    <- file.path(OUT_DIR, "univariate.tsv")
FN_SUM    <- file.path(OUT_DIR, "summary.txt")

# ======================
# helpers
# ======================
split_num <- function(x){
  if (is.na(x) || x=="") return(numeric())
  suppressWarnings(as.numeric(str_split(x, "[,;]")[[1]]))
}
sum_counts <- function(x){ v <- split_num(x); v[is.na(v)] <- 0; sum(v) }

# safe compact coordinates (optional, handy in outputs)
coord_string <- function(df_row){
  have <- function(nm) !is.null(df_row[[nm]]) && !is.na(df_row[[nm]])
  ev <- df_row$event; chr <- df_row$chr; strand <- df_row$strand
  if (ev=="SE"  && have("ExonStart_0base") && have("ExonEnd"))
    return(sprintf("%s:%s %d-%d", chr,strand, df_row$ExonStart_0base, df_row$ExonEnd))
  if (ev=="MXE" && all(c("X1stExonStart_0base","X1stExonEnd","X2ndExonStart_0base","X2ndExonEnd") %in% names(df_row)))
    return(sprintf("%s:%s E1:%d-%d vs E2:%d-%d", chr,strand,
                   df_row$X1stExonStart_0base, df_row$X1stExonEnd,
                   df_row$X2ndExonStart_0base, df_row$X2ndExonEnd))
  if (ev %in% c("A5SS","A3SS") && all(c("longExonStart_0base","longExonEnd","shortExonStart_0base","shortExonEnd") %in% names(df_row)))
    return(sprintf("%s:%s long:%d-%d; short:%d-%d", chr,strand,
                   df_row$longExonStart_0base, df_row$longExonEnd,
                   df_row$shortExonStart_0base, df_row$shortExonEnd))
  if (ev=="RI" && have("riExonStart_0base") && have("riExonEnd"))
    return(sprintf("%s:%s intron:%d-%d", chr,strand, df_row$riExonStart_0base, df_row$riExonEnd))
  NA_character_
}

read_one <- function(path, ev){
  if (!file.exists(path)) return(NULL)
  x <- suppressMessages(read_tsv(path, show_col_types = FALSE)) %>% mutate(event = ev)
  
  # base fields always kept
  base_cols <- c("event","geneSymbol","chr","strand","FDR","IncLevelDifference",
                 "IncLevel1","IncLevel2","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  # event-specific coordinates (optional)
  ev_cols <- switch(ev,
                    SE   = c("ExonStart_0base","ExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"),
                    MXE  = c("X1stExonStart_0base","X1stExonEnd","X2ndExonStart_0base","X2ndExonEnd",
                             "upstreamES","upstreamEE","downstreamES","downstreamEE"),
                    A5SS = c("longExonStart_0base","longExonEnd","shortExonStart_0base","shortExonEnd","flankingES","flankingEE"),
                    A3SS = c("longExonStart_0base","longExonEnd","shortExonStart_0base","shortExonEnd","flankingES","flankingEE"),
                    RI   = c("riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")
  )
  
  keep <- intersect(c(base_cols, ev_cols), names(x))
  x <- x %>% dplyr::select(all_of(keep))
  
  # quick filters to define “candidates”
  x <- x %>%
    rowwise() %>%
    dplyr::mutate(
      reads_HC    = sum_counts(IJC_SAMPLE_1) + sum_counts(SJC_SAMPLE_1),
      reads_MDD   = sum_counts(IJC_SAMPLE_2) + sum_counts(SJC_SAMPLE_2),
      reads_total = reads_HC + reads_MDD,
      coord       = coord_string(cur_data())
    ) %>% ungroup() %>%
    dplyr::filter(!is.na(FDR), !is.na(IncLevelDifference),
           FDR <= FDR_MAX,
           abs(IncLevelDifference) >= DPSI_MIN,
           reads_total >= READS_MIN)
  
  if (nrow(x) == 0) return(NULL)
  
  # build per-sample PSI vector for each row: c(IncLevel1 samples, IncLevel2 samples)
  n1 <- length(split_num(x$IncLevel1[1] %||% NA_character_))
  n2 <- length(split_num(x$IncLevel2[1] %||% NA_character_))
  if (n1 == 0 || n2 == 0) return(NULL)
  
  psi_mat <- matrix(NA_real_, nrow = nrow(x), ncol = n1 + n2)
  for (i in seq_len(nrow(x))) {
    v1 <- split_num(x$IncLevel1[i]); v2 <- split_num(x$IncLevel2[i])
    # make lengths consistent per row
    if (length(v1) != n1 || length(v2) != n2) next
    psi_mat[i, ] <- c(v1, v2)
  }
  
  # feature metadata
  feat <- x %>% transmute(
    feature_id = paste(event, geneSymbol, row_number(), sep="|"),
    event, geneSymbol, chr, strand, FDR, IncLevelDifference, reads_total, coord
  )
  
  list(psi = psi_mat, meta = feat, n1 = n1, n2 = n2)
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

# ======================
# LOAD & BUILD MATRIX
# ======================
lst <- imap(paths, read_one) %>% compact()
stopifnot(length(lst) > 0)

# ensure consistent group sizes across event types
n1s <- map_int(lst, "n1"); n2s <- map_int(lst, "n2")
if (length(unique(n1s)) != 1 || length(unique(n2s)) != 1)
  stop("Group sizes differ across files; make sure all rMATS runs used same samples/order.")

n1 <- unique(n1s); n2 <- unique(n2s)
y  <- c(rep(0L, n1), rep(1L, n2))  # 0=HC (IncLevel1), 1=MDD (IncLevel2)

# stack features
PSI <- do.call(rbind, lapply(lst, `[[`, "psi"))
META <- bind_rows(lapply(lst, `[[`, "meta"))
stopifnot(nrow(PSI) == nrow(META))

# drop features with too many NAs or near-zero variance
na_prop <- rowMeans(is.na(PSI))
keep_rows <- which(na_prop <= MAX_NA_PROP & apply(PSI, 1, function(v) sd(v, na.rm=TRUE) > 1e-6))
PSI <- PSI[keep_rows, , drop=FALSE]
META <- META[keep_rows, , drop=FALSE]
stopifnot(nrow(PSI) > 0)

# transpose to samples x features for glmnet
X <- t(PSI); colnames(X) <- META$feature_id; rownames(X) <- paste0("S", seq_len(nrow(X)))
write_tsv(as_tibble(X, rownames = "sample"), FN_PSI)
write_tsv(tibble(sample = rownames(X), label = y), FN_META)

# ======================
# UNIVARIATE SCREEN
# ======================
uni_p <- apply(X, 2, function(f){
  df <- data.frame(y=y, f=f)
  if (all(is.na(f))) return(1)
  # glm drops rows with NA automatically
  fit <- try(suppressWarnings(glm(y ~ f, data=df, family="binomial")), silent=TRUE)
  if (inherits(fit,"try-error")) return(1)
  co <- coef(summary(fit))
  if (nrow(co) < 2) return(1)
  co[2,4]
})
UNI <- tibble(feature_id = names(uni_p), p = uni_p) %>%
  arrange(p) %>% mutate(rank = row_number())
write_tsv(UNI, FN_UNI)

top_ids <- UNI$feature_id[seq_len(min(TOP_UNI, sum(is.finite(UNI$p))))]
X_top <- X[, top_ids, drop=FALSE]
META_top <- META %>% dplyr::filter(feature_id %in% top_ids)

# ======================
# MULTIVARIATE: LASSO (fallback to elastic-net if empty)
# ======================
fit_and_eval <- function(Xm, y, alpha=1){
  cv <- cv.glmnet(Xm, y, family="binomial", alpha=alpha, nfolds=N_FOLDS, standardize=TRUE)
  # prefer sparse s=lambda.1se, but try lambda.min if 0 selected
  s_use <- "lambda.1se"
  nz <- which(coef(cv, s=s_use)[-1,1] != 0)
  if (length(nz) == 0) {
    s_use <- "lambda.min"
    nz <- which(coef(cv, s=s_use)[-1,1] != 0)
  }
  pred <- as.numeric(predict(cv, newx = Xm, s = s_use, type = "response"))
  auc  <- as.numeric(pROC::roc(y, pred)$auc)
  list(cv=cv, s=s_use, nz=nz, pred=pred, auc=auc)
}

res <- fit_and_eval(X_top, y, alpha = USE_ALPHA)
if (length(res$nz) == 0 && USE_ALPHA == 1) {
  # fallback: elastic-net
  res <- fit_and_eval(X_top, y, alpha = 0.5)
}

# ======================
# REPORT
# ======================
sel <- tibble()
if (length(res$nz) > 0) {
  sel_idx <- res$nz
  sel_ids <- colnames(X_top)[sel_idx]
  coefs   <- as.numeric(coef(res$cv, s=res$s)[sel_idx+1,1])
  sel <- META_top %>%
    dplyr::mutate(feature_id = factor(feature_id, levels = colnames(X_top))) %>%
    dplyr::filter(feature_id %in% sel_ids) %>%
    dplyr::arrange(match(feature_id, sel_ids)) %>%
    dplyr::mutate(coef = coefs) %>%
    dplyr::select(feature_id, coef, event, geneSymbol, chr, strand, coord, FDR, IncLevelDifference, reads_total)
  write_tsv(sel, FN_PANEL)
}

sink(FN_SUM)
cat("Samples:", nrow(X), " | Features before uni:", ncol(X), " | After uni:", ncol(X_top), "\n")
cat("Selected:", nrow(sel), " | In-sample AUC:", round(res$auc, 3), " | alpha used:", ifelse(length(res$nz)>0, ifelse(res$cv$alpha==1,"1 (LASSO)","0.5 (EN)"), USE_ALPHA), "\n")
sink()

cat("Done.\n")
cat("Selected:", nrow(sel), " | In-sample AUC:", round(res$auc, 3), "\n")
cat("Files:\n  - ", FN_PSI, "\n  - ", FN_META, "\n  - ", FN_UNI, "\n  - ", FN_PANEL, "\n  - ", FN_SUM, "\n", sep="")
