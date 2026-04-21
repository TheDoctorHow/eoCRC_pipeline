##############################################################################
# step6_network_analysis.R
# 
# eoCRC Network Analysis Pipeline
# Builds SpiecEasi co-occurrence networks, runs community detection,
# network diffusion, dual-module scoring, and patient stratification.
#
# Inputs:  species_filtered.rds, all_samples_metadata.rds, tse_pathways_full.rds
# Outputs: results/networks/ (all network objects, scores, summaries)
#
# Runtime: ~30-40 minutes (SpiecEasi is the bottleneck)
##############################################################################

# ── Portable directory setup ─────────────────────────────────────
if (!exists("BASE_DIR")) {
  BASE_DIR <- normalizePath("~/eoCRC_analysis")
}
NETWORK_DIR <- file.path(BASE_DIR, "results", "networks")
dir.create(NETWORK_DIR, recursive = TRUE, showWarnings = FALSE)

cat("══════════════════════════════════════════\n")
cat("  STEP 6: NETWORK ANALYSIS\n")
cat("══════════════════════════════════════════\n\n")
cat("Base directory:", BASE_DIR, "\n")
cat("Output directory:", NETWORK_DIR, "\n\n")

# ── Load packages ────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(SpiecEasi)
  library(igraph)
  library(pROC)
  library(vegan)
  library(ppcor)
})

# ── Load source data ─────────────────────────────────────────────
cat("Loading data...\n")
sp <- readRDS(file.path(BASE_DIR, "species_filtered.rds"))
meta_all <- readRDS(file.path(BASE_DIR, "all_samples_metadata.rds"))
rownames(meta_all) <- meta_all$sample_id
meta_all$four_group <- paste0(
  ifelse(meta_all$group == "CRC", "CRC", "Ctrl"),
  "_",
  ifelse(meta_all$age < 50, "young", "older")
)
cat("Species matrix:", dim(sp), "\n")
cat("Groups:", table(meta_all$four_group), "\n\n")

# ── Species filter (>20% prevalence + forced oral pathobionts) ───
prev_all <- colMeans(sp > 0)
oral_must_include <- grep("Peptostreptococcus stomatis|Dialister pneumosintes",
                          colnames(sp), value = TRUE)
keep_cols <- unique(c(colnames(sp)[prev_all > 0.20], oral_must_include))
sp_filt <- sp[, keep_cols]
short_names <- gsub("species:", "", colnames(sp_filt))
cat("Filtered species:", ncol(sp_filt), "\n\n")

# ── Define species lists ─────────────────────────────────────────
oral_taxa <- c("Parvimonas micra", "Gemella morbillorum",
               "Peptostreptococcus stomatis", "Dialister pneumosintes")

crc_associated <- c("Bacteroides fragilis", "Fusobacterium nucleatum",
                    "[Clostridium] symbiosum", "Enterocloster clostridioformis",
                    "Eisenbergiella tayi", "Hungatella hathewayi",
                    "Porphyromonas asaccharolytica", "[Clostridium] innocuum")

crc_taxa_broad <- c("Bacteroides fragilis", "Fusobacterium nucleatum",
                    "[Clostridium] symbiosum", "Enterocloster clostridioformis",
                    "Eisenbergiella tayi", "Hungatella hathewayi",
                    "[Clostridium] innocuum")

butyrate_producers <- c("Faecalibacterium prausnitzii", "Roseburia intestinalis",
                        "Roseburia inulinivorans", "Anaerostipes hadrus",
                        "Coprococcus eutactus")

oral_origin <- grep("Parvimonas|Gemella|Peptostreptococcus stomatis|Dialister pneumosintes|Fusobacterium nucleatum|Solobacterium|Mogibacterium|Eubacterium sulci|Streptococcus gordonii|Streptococcus sanguinis|Actinomyces|Veillonella",
                    short_names, value = TRUE)

saveRDS(list(oral_taxa = oral_taxa, oral_origin = oral_origin,
             crc_associated = crc_associated, crc_taxa_broad = crc_taxa_broad,
             butyrate_producers = butyrate_producers, short_names = short_names),
        file.path(NETWORK_DIR, "species_annotations.rds"))

# ══════════════════════════════════════════════════════════════════
# 1. BUILD SpiecEasi NETWORKS
# ══════════════════════════════════════════════════════════════════

cat("══ Building SpiecEasi networks ══\n\n")

build_network <- function(group_name, sp_filt, meta_all, short_names, NETWORK_DIR) {
  idx <- meta_all$four_group == group_name
  cat("  ", group_name, "(n=", sum(idx), ")...")
  se <- spiec.easi(sp_filt[idx, ], method = "mb",
                    lambda.min.ratio = 1e-2, nlambda = 20,
                    pulsar.params = list(rep.num = 50, seed = 42))
  adj <- as.matrix(getRefit(se))
  colnames(adj) <- rownames(adj) <- short_names
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  cat(" edges:", ecount(g), "\n")
  saveRDS(se, file.path(NETWORK_DIR, paste0("spieceasi_", group_name, ".rds")))
  return(g)
}

g_eoCRC <- build_network("CRC_young", sp_filt, meta_all, short_names, NETWORK_DIR)
g_yCtrl <- build_network("Ctrl_young", sp_filt, meta_all, short_names, NETWORK_DIR)
g_loCRC <- build_network("CRC_older", sp_filt, meta_all, short_names, NETWORK_DIR)
g_oCtrl <- build_network("Ctrl_older", sp_filt, meta_all, short_names, NETWORK_DIR)

graphs <- list(eoCRC = g_eoCRC, yCtrl = g_yCtrl, loCRC = g_loCRC, oCtrl = g_oCtrl)
saveRDS(graphs, file.path(NETWORK_DIR, "igraph_objects.rds"))

# ══════════════════════════════════════════════════════════════════
# 2. NETWORK DIFFUSION FUNCTION
# ══════════════════════════════════════════════════════════════════

diffuse_from <- function(g, seed_species, restart_prob = 0.15, n_iter = 100) {
  n <- vcount(g)
  node_names <- V(g)$name
  seed_idx <- which(node_names %in% seed_species)
  if (length(seed_idx) == 0) return(NULL)
  p <- rep(0, n)
  p[seed_idx] <- 1 / length(seed_idx)
  p0 <- p
  A <- as.matrix(as_adjacency_matrix(g))
  col_sums <- colSums(A)
  col_sums[col_sums == 0] <- 1
  W <- t(t(A) / col_sums)
  for (i in 1:n_iter) {
    p <- (1 - restart_prob) * (W %*% p) + restart_prob * p0
  }
  result <- data.frame(species = node_names, influence = as.numeric(p),
                       stringsAsFactors = FALSE)
  result[order(-result$influence), ]
}

compute_ratio <- function(g, seeds, crc_list, bp_list) {
  d <- diffuse_from(g, seeds)
  d <- d[!d$species %in% seeds, ]
  crc_inf <- sum(d$influence[d$species %in% crc_list])
  bp_inf <- sum(d$influence[d$species %in% bp_list])
  list(ratio = crc_inf / max(bp_inf, 1e-6),
       crc_influence = crc_inf, bp_influence = bp_inf, diffusion = d)
}

# ══════════════════════════════════════════════════════════════════
# 3. COMPUTE DIFFUSION RATIOS
# ══════════════════════════════════════════════════════════════════

cat("\n══ Network diffusion analysis ══\n\n")

diff_results <- list()
for (nm in names(graphs)) {
  r <- compute_ratio(graphs[[nm]], oral_taxa, crc_taxa_broad, butyrate_producers)
  diff_results[[nm]] <- r
  cat(nm, ": CRC:BP ratio =", round(r$ratio, 2), "\n")
}

# ── Permutation test (eoCRC) ─────────────────────────────────────
cat("\nRunning permutation test (1000 permutations)...\n")
set.seed(42)
n_perm <- 1000
non_target <- setdiff(V(g_eoCRC)$name, c(oral_taxa, crc_taxa_broad, butyrate_producers))
perm_ratios <- numeric(n_perm)
for (i in 1:n_perm) {
  rand_seeds <- sample(non_target, 4)
  r <- compute_ratio(g_eoCRC, rand_seeds, crc_taxa_broad, butyrate_producers)
  perm_ratios[i] <- r$ratio
}
perm_p <- mean(perm_ratios >= diff_results$eoCRC$ratio)
cat("Observed ratio:", round(diff_results$eoCRC$ratio, 2),
    "| Permutation p:", perm_p,
    "| Mean random:", round(mean(perm_ratios), 2), "\n")

saveRDS(list(
  ratios = sapply(diff_results, function(x) x$ratio),
  permutation_p = perm_p,
  permutation_distribution = perm_ratios,
  diffusions = lapply(diff_results, function(x) x$diffusion)
), file.path(NETWORK_DIR, "diffusion_results.rds"))

# ══════════════════════════════════════════════════════════════════
# 4. COMMUNITY DETECTION
# ══════════════════════════════════════════════════════════════════

cat("\n══ Community detection (100 Louvain iterations) ══\n\n")

set.seed(42)
n_runs <- 100
stability_matrix <- matrix(0, nrow = n_runs, ncol = 4,
                           dimnames = list(NULL, names(graphs)))

for (i in 1:n_runs) {
  for (j in seq_along(graphs)) {
    comm <- cluster_louvain(graphs[[j]])
    mem <- membership(comm)
    count <- 0
    for (bug in oral_taxa) {
      comm_members <- names(mem[mem == mem[bug]])
      count <- count + length(intersect(comm_members, crc_associated))
    }
    stability_matrix[i, j] <- count
  }
}

cat("Mean co-clusterings:\n")
print(round(colMeans(stability_matrix), 1))
cat("Proportion > 0:\n")
print(round(colMeans(stability_matrix > 0), 2))

saveRDS(stability_matrix, file.path(NETWORK_DIR, "community_stability_matrix.rds"))

# ══════════════════════════════════════════════════════════════════
# 5. DUAL-MODULE SCORES
# ══════════════════════════════════════════════════════════════════

cat("\n══ Dual-module scores ══\n\n")

oral_crc_cols <- sapply(c("Parvimonas micra", "Peptostreptococcus stomatis",
                          "Gemella morbillorum", "Bacteroides fragilis",
                          "Fusobacterium nucleatum", "Hungatella hathewayi"),
                        function(s) grep(s, colnames(sp), value = TRUE)[1])
butyrate_cols <- sapply(butyrate_producers,
                        function(s) grep(s, colnames(sp), value = TRUE)[1])
butyrate_cols <- butyrate_cols[!is.na(butyrate_cols)]

sp_pseudo <- sp + 0.01
sp_clr <- t(apply(sp_pseudo, 1, function(x) log(x) - mean(log(x))))

oralcrc_score <- rowMeans(sp_clr[, oral_crc_cols])
bp_score <- rowMeans(sp_clr[, butyrate_cols])
dual_ratio <- oralcrc_score - bp_score

for (g in c("CRC_young", "Ctrl_young", "CRC_older", "Ctrl_older")) {
  idx <- meta_all$four_group == g
  cat(g, ": dual ratio mean =", round(mean(dual_ratio[idx]), 3), "\n")
}

# Anticorrelation
cat("\nModule anticorrelation (Spearman rho):\n")
for (g in c("CRC_young", "Ctrl_young", "CRC_older", "Ctrl_older")) {
  idx <- meta_all$four_group == g
  ct <- cor.test(oralcrc_score[idx], bp_score[idx], method = "spearman")
  cat("  ", g, ": rho =", round(ct$estimate, 3), " p =", signif(ct$p.value, 3), "\n")
}

# Classification
cat("\nDual ratio AUC (all CRC vs controls):",
    round(as.numeric(auc(roc(as.numeric(meta_all$group == "CRC"),
                             dual_ratio, quiet = TRUE))), 3), "\n")

module_scores <- data.frame(
  sample_id = meta_all$sample_id,
  four_group = meta_all$four_group,
  oralcrc_score = oralcrc_score,
  bp_score = bp_score,
  dual_ratio = dual_ratio
)
write.csv(module_scores, file.path(NETWORK_DIR, "module_scores.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════
# 6. WITHIN-STUDY REPRODUCIBILITY
# ══════════════════════════════════════════════════════════════════

cat("\n══ Within-study co-occurrence reproducibility ══\n\n")

eoCRC_idx <- meta_all$four_group == "CRC_young"

pairs <- list(
  c("Gemella morbillorum", "Peptostreptococcus stomatis"),
  c("Parvimonas micra", "Fusobacterium nucleatum"),
  c("Peptostreptococcus stomatis", "Bacteroides fragilis"),
  c("Parvimonas micra", "Bacteroides fragilis"),
  c("Peptostreptococcus stomatis", "Fusobacterium nucleatum")
)

reproducibility_results <- list()
for (p in pairs) {
  col1 <- grep(p[1], colnames(sp), value = TRUE)[1]
  col2 <- grep(p[2], colnames(sp), value = TRUE)[1]
  pair_name <- paste(p[1], "<->", p[2])
  
  study_results <- list()
  for (study in unique(meta_all$study_name[eoCRC_idx])) {
    idx <- eoCRC_idx & meta_all$study_name == study
    if (sum(idx) >= 6) {
      ct <- cor.test(sp[idx, col1], sp[idx, col2], method = "spearman")
      study_results[[study]] <- list(rho = ct$estimate, p = ct$p.value, n = sum(idx))
    }
  }
  
  signs <- sapply(study_results, function(x) sign(x$rho))
  cat(pair_name, ":", sum(signs > 0), "positive /", sum(signs < 0), "negative of",
      length(signs), "studies\n")
  reproducibility_results[[pair_name]] <- study_results
}

# Partial correlations controlling for study
cat("\nPartial correlations (controlling for study):\n")
eo_data <- data.frame(
  ps = sp[eoCRC_idx, grep("Peptostreptococcus stomatis", colnames(sp), value = TRUE)[1]],
  bf = sp[eoCRC_idx, grep("Bacteroides fragilis", colnames(sp), value = TRUE)[1]],
  fn = sp[eoCRC_idx, grep("Fusobacterium nucleatum", colnames(sp), value = TRUE)[1]],
  pm = sp[eoCRC_idx, grep("Parvimonas micra", colnames(sp), value = TRUE)[1]],
  study = as.numeric(as.factor(meta_all$study_name[eoCRC_idx]))
)

partial_pairs <- list(
  c("ps", "bf", "P.stomatis <-> B.fragilis"),
  c("pm", "fn", "P.micra <-> F.nucleatum"),
  c("ps", "fn", "P.stomatis <-> F.nucleatum")
)

for (pp in partial_pairs) {
  pc <- pcor.test(eo_data[[pp[1]]], eo_data[[pp[2]]], eo_data$study, method = "spearman")
  cat("  ", pp[3], ": partial rho =", round(pc$estimate, 3),
      " p =", signif(pc$p.value, 3), "\n")
}

saveRDS(reproducibility_results, file.path(NETWORK_DIR, "reproducibility_results.rds"))

# ══════════════════════════════════════════════════════════════════
# 7. PATIENT STRATIFICATION
# ══════════════════════════════════════════════════════════════════

cat("\n══ Patient stratification ══\n\n")

oral_cols_sp <- sapply(c("Parvimonas micra", "Gemella morbillorum",
                          "Peptostreptococcus stomatis"),
                        function(s) grep(s, colnames(sp), value = TRUE)[1])
crc_cols_sp <- sapply(c("Bacteroides fragilis", "Hungatella hathewayi",
                         "Fusobacterium nucleatum"),
                       function(s) grep(s, colnames(sp), value = TRUE)[1])

oral_count <- rowSums(sp[, oral_cols_sp] > 0)
crc_count <- rowSums(sp[, crc_cols_sp] > 0)

eo_integrated <- oral_count[eoCRC_idx] >= 2 & crc_count[eoCRC_idx] >= 2
lo_integrated <- oral_count[meta_all$four_group == "CRC_older"] >= 2 &
                 crc_count[meta_all$four_group == "CRC_older"] >= 2

cat("Integrated patients:\n")
cat("  eoCRC:", sum(eo_integrated), "/", sum(eoCRC_idx),
    "(", round(100 * mean(eo_integrated), 1), "%)\n")
cat("  loCRC:", sum(lo_integrated), "/", sum(meta_all$four_group == "CRC_older"),
    "(", round(100 * mean(lo_integrated), 1), "%)\n")

# Dual ratio in integrated vs non-integrated
eo_dual <- dual_ratio[eoCRC_idx]
cat("\nDual ratio — integrated vs non-integrated (eoCRC):\n")
cat("  Integrated:", round(mean(eo_dual[eo_integrated]), 3), "\n")
cat("  Non-integrated:", round(mean(eo_dual[!eo_integrated]), 3), "\n")
cat("  Wilcoxon p:", signif(wilcox.test(eo_dual[eo_integrated],
                                        eo_dual[!eo_integrated])$p.value, 3), "\n")

# ══════════════════════════════════════════════════════════════════
# 8. NETWORK SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════

network_summary <- data.frame(
  Group = c("eoCRC", "Young Ctrl", "loCRC", "Older Ctrl"),
  n = c(78, 121, 562, 521),
  total_edges = sapply(graphs, ecount),
  density = round(sapply(graphs, edge_density), 4),
  transitivity = round(sapply(graphs, function(g) transitivity(g, type = "global")), 4),
  diffusion_ratio = round(sapply(diff_results, function(x) x$ratio), 2),
  mean_coclusterings = round(colMeans(stability_matrix), 1)
)
write.csv(network_summary, file.path(NETWORK_DIR, "network_summary.csv"), row.names = FALSE)
cat("\n")
print(network_summary)

cat("\n══ Step 6 complete. Results saved to", NETWORK_DIR, "══\n")
