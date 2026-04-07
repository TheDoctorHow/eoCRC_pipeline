#!/usr/bin/env Rscript
# Figure 6 — Pathway heatmap
# Fixes: (1) clean cohort legend labels, (2) fixed pathway label formatting,
#        (3) q-values from MaAsLin2 v3, (4) 300 dpi PNG + cairo_pdf

suppressPackageStartupMessages({
  library(dplyr); library(compositions); library(pheatmap)
})

if (basename(getwd()) == "scripts") setwd("..")
dir.create("figures", showWarnings = FALSE)

if (!file.exists("all_samples_metadata.rds") ||
    !file.exists("path_filtered.rds"))
  stop("Run scripts/step1_extract_filter.R first to generate input data.")
if (!file.exists("maaslin2_results_v2/summaries/all_associations.csv"))
  stop("Run scripts/step2b_maaslin2_corrected.R first to generate MaAsLin2 v2 results.")
if (!file.exists("maaslin2_results_v3/summaries/all_associations.csv"))
  stop("Run scripts/step2c_maaslin2_v3.R first to generate MaAsLin2 v3 results.")

# ── Cohort colour palette (raw study names → colours) ─────────────────────────
COHORT_COLS_RAW <- c(
  FengQ_2015      = "#E67E22", GuptaA_2019     = "#16A085",
  HanniganGD_2017 = "#8E44AD", ThomasAM_2019_c = "#2980B9",
  VogtmannE_2016  = "#C0392B", WirbelJ_2018    = "#27AE60",
  YachidaS_2019   = "#D35400", YuJ_2015        = "#2C3E50",
  ZellerG_2014    = "#E91E63"
)

# Clean cohort labels for the legend
COHORT_CLEAN <- c(
  FengQ_2015      = "Feng 2015",    GuptaA_2019     = "Gupta 2019",
  HanniganGD_2017 = "Hannigan 2017", ThomasAM_2019_c = "Thomas 2019",
  VogtmannE_2016  = "Vogtmann 2016", WirbelJ_2018    = "Wirbel 2018",
  YachidaS_2019   = "Yachida 2019", YuJ_2015        = "Yu 2015",
  ZellerG_2014    = "Zeller 2014"
)

COL_eo  <- "#E74C3C"
COL_yc  <- "#27AE60"

# ── Clean pathway label lookup ─────────────────────────────────────────────────
# Keys = matrix column names (as in path_filtered.rds); values = display labels
PW_LABELS <- c(
  # User-specified fixes
  "FASYN-ELONG-PWY: fatty acid elongation -- saturated"               = "fatty acid elongation (saturated)",
  "PWY-6282: palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)"  = "palmitoleate biosynthesis I (from 5Z-dodec-5-enoate)",
  "PWY0-862: (5Z)-dodec-5-enoate biosynthesis"                        = "5Z-dodec-5-enoate biosynthesis",
  "PWY-6519: 8-amino-7-oxononanoate biosynthesis I"                   = "8-amino-7-oxononanoate biosynthesis I",
  "PWY-7237: myo-, chiro- and scillo-inositol degradation"            = "myo-, chiro-, and scillo-inositol degradation",
  # Other top-20 pathways (clean names from matrix, but kept for completeness)
  "THISYNARA-PWY: superpathway of thiamin diphosphate biosynthesis III (eukaryotes)" =
    "thiamin diphosphate biosynthesis III (superpathway)",
  # Replacement pathway (backfill for ppGpp biosynthesis)
  "P163-PWY: L-lysine fermentation to acetate and butanoate"                        =
    "L-lysine fermentation to acetate and butanoate",
  "PWYG-321: mycolate biosynthesis"                                   = "mycolate biosynthesis",
  "PWY-7234: inosine-5'-phosphate biosynthesis III"                   = "inosine-5'-phosphate biosynthesis III",
  "PWY-5695: urate biosynthesis/inosine 5'-phosphate degradation"     = "urate biosynthesis / inosine-5'-P degradation",
  "PWY-5989: stearate biosynthesis II (bacteria and plants)"          = "stearate biosynthesis II (bacteria and plants)"
)

# Generic label extractor: strip PWY-ID prefix, then apply lookup / fallback
clean_pw_label <- function(col_name) {
  # First check exact lookup
  if (col_name %in% names(PW_LABELS)) return(PW_LABELS[[col_name]])
  # Otherwise: take description after ": "
  desc <- sub("^[^:]+:\\s*", "", col_name)
  desc
}

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading data...\n")
meta <- readRDS("all_samples_metadata.rds") %>%
  mutate(grp = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))
pa_mat  <- readRDS("path_filtered.rds")
maas_v2 <- read.csv("maaslin2_results_v2/summaries/all_associations.csv",
                    stringsAsFactors = FALSE)
maas_v3 <- read.csv("maaslin2_results_v3/summaries/all_associations.csv",
                    stringsAsFactors = FALSE)

pseudo  <- 1e-6

# ── Primary samples (eoCRC + young ctrl) ──────────────────────────────────────
idx_pri  <- meta$grp %in% c("eoCRC", "young_ctrl")
pa_pri   <- pa_mat[idx_pri, ]
pa_clr_p <- as.matrix(clr(pa_pri + pseudo))
colnames(pa_clr_p) <- colnames(pa_mat)   # preserve column names after CLR

# Helper: convert MaAsLin2 feature name (dots) → PWY ID (hyphens)
get_pwy_id <- function(feat) {
  id_part <- sub("\\.\\..*", "", feat)   # before first double-dot
  gsub("\\.", "-", id_part)
}

# All v3 pathway results sorted by q-value (used for both filtering & backfill)
v3_pw_all <- maas_v3 %>%
  filter(comparison == "eoCRC_vs_YoungCtrl", feature_set == "pathways") %>%
  arrange(qval, desc(abs(coef))) %>%
  mutate(pwy_id = get_pwy_id(feature))

# Match a MaAsLin2 feature to a matrix column by PWY-ID prefix
pa_colnames <- colnames(pa_clr_p)
match_col <- function(feat) {
  id  <- get_pwy_id(feat)
  hit <- pa_colnames[startsWith(pa_colnames, paste0(id, ":"))]
  if (length(hit) > 0) hit[1] else NA_character_
}

# ── Build initial set: v2 top-20 with v3 q-values attached ───────────────────
top20_v2 <- maas_v2 %>%
  filter(comparison == "eoCRC_vs_YoungCtrl", feature_set == "pathways") %>%
  arrange(qval, desc(abs(coef))) %>%
  slice_head(n = 20) %>%
  mutate(mat_col = sapply(feature, match_col)) %>%
  filter(!is.na(mat_col))

v3_pw <- v3_pw_all %>% select(feature, coef_v3 = coef, qval_v3 = qval)
top20_v2 <- left_join(top20_v2, v3_pw, by = "feature")

# ── Filter: remove pathways with v3 q ≥ 0.25; backfill from v3 ───────────────
removed <- top20_v2 %>% filter(!is.na(qval_v3) & qval_v3 >= 0.25)
if (nrow(removed) > 0) {
  cat(sprintf("Removing %d pathway(s) with v3 q >= 0.25:\n", nrow(removed)))
  for (i in seq_len(nrow(removed)))
    cat(sprintf("  %s (v3 q=%.3f)\n", removed$feature[i], removed$qval_v3[i]))
}

top20_v2 <- top20_v2 %>% filter(is.na(qval_v3) | qval_v3 < 0.25)

# Backfill with top v3 candidates not already present, q < 0.25, in matrix
already_ids <- get_pwy_id(top20_v2$feature)
backfill <- v3_pw_all %>%
  filter(qval < 0.25, !is.na(pwy_id), !(pwy_id %in% already_ids)) %>%
  mutate(mat_col = sapply(feature, match_col)) %>%
  filter(!is.na(mat_col)) %>%
  slice_head(n = max(0, 20 - nrow(top20_v2)))

if (nrow(backfill) > 0) {
  cat(sprintf("Adding %d replacement pathway(s) from v3:\n", nrow(backfill)))
  for (i in seq_len(nrow(backfill)))
    cat(sprintf("  %s (v3 q=%.4f)\n", backfill$feature[i], backfill$qval[i]))
  backfill <- backfill %>%
    rename(coef_v3 = coef, qval_v3 = qval) %>%
    mutate(coef = coef_v3, qval = qval_v3)   # use v3 values in coef/qval slots
  top20_v2 <- bind_rows(top20_v2, backfill)
}

cat(sprintf("Final pathway count: %d\n", nrow(top20_v2)))

# ── Build heatmap matrix ───────────────────────────────────────────────────────
heat_mat <- pa_clr_p[, top20_v2$mat_col, drop = FALSE]   # samples × pathways

# Build row labels: "clean name [↑/↓, q=X.XXX]" using v3 values
make_row_label <- function(mat_col, coef_v3, qval_v3) {
  desc <- clean_pw_label(mat_col)
  dir  <- if (!is.na(coef_v3) && coef_v3 > 0) "\u2191" else "\u2193"
  q_str <- if (is.na(qval_v3)) "q=NA"
            else if (qval_v3 < 0.001) sprintf("q=%s", formatC(qval_v3, format = "e", digits = 1))
            else sprintf("q=%.3f", qval_v3)
  sprintf("%s [%s, %s]", desc, dir, q_str)
}

row_labels <- mapply(make_row_label,
                     top20_v2$mat_col, top20_v2$coef_v3, top20_v2$qval_v3)

# ── Column annotation: Group + Cohort (clean labels) ──────────────────────────
study_raw   <- meta$study_name[idx_pri]
study_clean <- COHORT_CLEAN[study_raw]   # recode to clean names

col_ann <- data.frame(
  Group  = ifelse(meta$grp[idx_pri] == "eoCRC", "eoCRC", "Young ctrl"),
  Cohort = study_clean,
  row.names = rownames(heat_mat)
)

# Colours keyed by clean cohort names
cohorts_present  <- unique(study_clean)
cohort_cols_clean <- setNames(
  COHORT_COLS_RAW[names(COHORT_CLEAN)[match(cohorts_present, COHORT_CLEAN)]],
  cohorts_present
)

ann_colors <- list(
  Group  = c(eoCRC = COL_eo, "Young ctrl" = COL_yc),
  Cohort = cohort_cols_clean
)

# ── Draw heatmap ──────────────────────────────────────────────────────────────
draw_heatmap <- function(fontsize_r = 7) {
  pheatmap(
    t(heat_mat),                      # rows = pathways, cols = samples
    color             = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    annotation_col    = col_ann,
    annotation_colors = ann_colors,
    show_colnames     = FALSE,
    fontsize_row      = fontsize_r,
    labels_row        = row_labels,
    main              = "Top 20 differentially abundant pathways \u2014 eoCRC vs Young Controls\nMaAsLin2 v3 CLR | direction and q-value annotated",
    border_color      = NA,
    scale             = "row"
  )
}

# PNG at 300 dpi
png("figures/fig6_pathway_heatmap.png", width = 3600, height = 3000, res = 300)
draw_heatmap(fontsize_r = 7)
invisible(dev.off())

# PDF via cairo_pdf for Unicode support
cairo_pdf("figures/fig6_pathway_heatmap.pdf", width = 12, height = 10)
draw_heatmap(fontsize_r = 7)
invisible(dev.off())

cat("Saved: figures/fig6_pathway_heatmap.png + .pdf\n")
cat("Done.\n")
