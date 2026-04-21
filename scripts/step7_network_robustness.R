##############################################################################
# step7_network_robustness.R
#
# Robustness and sensitivity analyses for the network diffusion findings.
# Must run after step6_network_analysis.R
#
# Inputs:  results/networks/ (from step6), species_filtered.rds, all_samples_metadata.rds
# Outputs: results/networks/ (additional robustness files)
#
# Runtime: ~3-4 hours (100-iteration matched comparison is the bottleneck)
##############################################################################

# ── Portable directory setup ─────────────────────────────────────
if (!exists("BASE_DIR")) {
  BASE_DIR <- normalizePath("~/eoCRC_analysis")
}
NETWORK_DIR <- file.path(BASE_DIR, "results", "networks")

cat("══════════════════════════════════════════\n")
cat("  STEP 7: NETWORK ROBUSTNESS CHECKS\n")
cat("══════════════════════════════════════════\n\n")

# ── Load packages ────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(SpiecEasi)
  library(igraph)
  library(vegan)
})

# ── Load data and step6 outputs ──────────────────────────────────
sp <- readRDS(file.path(BASE_DIR, "species_filtered.rds"))
meta_all <- readRDS(file.path(BASE_DIR, "all_samples_metadata.rds"))
rownames(meta_all) <- meta_all$sample_id
meta_all$four_group <- paste0(
  ifelse(meta_all$group == "CRC", "CRC", "Ctrl"),
  "_",
  ifelse(meta_all$age < 50, "young", "older")
)

prev_all <- colMeans(sp > 0)
oral_must_include <- grep("Peptostreptococcus stomatis|Dialister pneumosintes",
                          colnames(sp), value = TRUE)
keep_cols <- unique(c(colnames(sp)[prev_all > 0.20], oral_must_include))
sp_filt <- sp[, keep_cols]
short_names <- gsub("species:", "", colnames(sp_filt))

annotations <- readRDS(file.path(NETWORK_DIR, "species_annotations.rds"))
oral_taxa <- annotations$oral_taxa
crc_taxa_broad <- annotations$crc_taxa_broad
butyrate_producers <- annotations$butyrate_producers
crc_associated <- annotations$crc_associated

graphs <- readRDS(file.path(NETWORK_DIR, "igraph_objects.rds"))
g_eoCRC <- graphs$eoCRC

# ── Diffusion function ───────────────────────────────────────────
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
  crc_inf / max(bp_inf, 1e-6)
}

# ══════════════════════════════════════════════════════════════════
# 1. LEAVE-ONE-STUDY-OUT DIFFUSION
# ══════════════════════════════════════════════════════════════════

cat("══ Leave-one-study-out diffusion ══\n\n")

eoCRC_idx <- meta_all$four_group == "CRC_young"
loso_ratios <- c()
loso_names <- c()

for (drop_study in unique(meta_all$study_name[eoCRC_idx])) {
  idx <- eoCRC_idx & meta_all$study_name != drop_study
  n_remain <- sum(idx)
  
  if (n_remain >= 30) {
    se_loo <- spiec.easi(sp_filt[idx, ], method = "mb",
                          lambda.min.ratio = 1e-2, nlambda = 20,
                          pulsar.params = list(rep.num = 50, seed = 42))
    adj_loo <- as.matrix(getRefit(se_loo))
    colnames(adj_loo) <- rownames(adj_loo) <- short_names
    g_loo <- graph_from_adjacency_matrix(adj_loo, mode = "undirected", diag = FALSE)
    
    ratio <- compute_ratio(g_loo, oral_taxa, crc_taxa_broad, butyrate_producers)
    loso_ratios <- c(loso_ratios, ratio)
    loso_names <- c(loso_names, drop_study)
    cat("  Drop", drop_study, "(n=", n_remain, "): ratio =", round(ratio, 2), "\n")
  }
}

# ══════════════════════════════════════════════════════════════════
# 2. COMMUNITY DETECTION — ALGORITHM ROBUSTNESS
# ══════════════════════════════════════════════════════════════════

cat("\n══ Community detection algorithm robustness ══\n\n")

algorithms <- list(
  louvain = function(g) cluster_louvain(g),
  walktrap = function(g) cluster_walktrap(g),
  infomap = function(g) cluster_infomap(g),
  label_prop = function(g) cluster_label_prop(g)
)

algo_results <- list()
for (group_name in c("eoCRC", "yCtrl")) {
  g <- graphs[[group_name]]
  algo_results[[group_name]] <- list()
  
  for (alg_name in names(algorithms)) {
    counts <- numeric(50)
    for (i in 1:50) {
      comm <- algorithms[[alg_name]](g)
      mem <- membership(comm)
      count <- 0
      for (bug in oral_taxa) {
        comm_members <- names(mem[mem == mem[bug]])
        count <- count + length(intersect(comm_members, crc_associated))
      }
      counts[i] <- count
    }
    algo_results[[group_name]][[alg_name]] <- counts
    cat(group_name, "-", alg_name, ": mean =", round(mean(counts), 1),
        " prop>0 =", round(mean(counts > 0), 2), "\n")
  }
}

# ══════════════════════════════════════════════════════════════════
# 3. SPECIES DEFINITION SENSITIVITY
# ══════════════════════════════════════════════════════════════════

cat("\n══ Species definition sensitivity ══\n\n")

crc_strict <- c("Bacteroides fragilis", "Fusobacterium nucleatum",
                "Hungatella hathewayi")
crc_broad <- crc_taxa_broad
bp_strict <- c("Faecalibacterium prausnitzii", "Roseburia intestinalis")
bp_broad <- butyrate_producers

sensitivity_results <- data.frame()
for (crc_label in c("strict_CRC", "broad_CRC")) {
  crc_set <- if (crc_label == "strict_CRC") crc_strict else crc_broad
  for (bp_label in c("strict_BP", "broad_BP")) {
    bp_set <- if (bp_label == "strict_BP") bp_strict else bp_broad
    
    r_eo <- compute_ratio(graphs$eoCRC, oral_taxa, crc_set, bp_set)
    r_ct <- compute_ratio(graphs$yCtrl, oral_taxa, crc_set, bp_set)
    
    cat(crc_label, "+", bp_label, ": eoCRC =", round(r_eo, 1),
        " yCtrl =", round(r_ct, 1), " fold =", round(r_eo / max(r_ct, 0.01), 1), "\n")
    
    sensitivity_results <- rbind(sensitivity_results, data.frame(
      crc_def = crc_label, bp_def = bp_label,
      eoCRC = r_eo, yCtrl = r_ct, fold_diff = r_eo / max(r_ct, 0.01)))
  }
}

# ══════════════════════════════════════════════════════════════════
# 4. RESTART PROBABILITY SENSITIVITY
# ══════════════════════════════════════════════════════════════════

cat("\n══ Restart probability sensitivity ══\n\n")

restart_results <- data.frame()
for (rp in c(0.05, 0.10, 0.15, 0.25, 0.50)) {
  ratios <- sapply(names(graphs), function(nm) {
    d <- diffuse_from(graphs[[nm]], oral_taxa, restart_prob = rp)
    d <- d[!d$species %in% oral_taxa, ]
    sum(d$influence[d$species %in% crc_taxa_broad]) /
      max(sum(d$influence[d$species %in% butyrate_producers]), 1e-6)
  })
  cat("restart =", rp, ":", paste(names(ratios), "=", round(ratios, 1), collapse = "  "), "\n")
  restart_results <- rbind(restart_results, data.frame(
    restart = rp, t(ratios)))
}

# ══════════════════════════════════════════════════════════════════
# 5. SAMPLE-SIZE-MATCHED loCRC (100 iterations)
# ══════════════════════════════════════════════════════════════════

cat("\n══ Sample-size-matched loCRC (100 iterations) ══\n")
cat("This will take 1-2 hours...\n\n")

set.seed(42)
n_iter <- 100
loCRC_matched_ratios <- numeric(n_iter)
loCRC_idx_all <- which(meta_all$four_group == "CRC_older")

for (b in 1:n_iter) {
  idx_sub <- sample(loCRC_idx_all, 78)
  se_sub <- spiec.easi(sp_filt[idx_sub, ], method = "mb",
                        lambda.min.ratio = 1e-2, nlambda = 20,
                        pulsar.params = list(rep.num = 50, seed = 42))
  adj_sub <- as.matrix(getRefit(se_sub))
  colnames(adj_sub) <- rownames(adj_sub) <- short_names
  g_sub <- graph_from_adjacency_matrix(adj_sub, mode = "undirected", diag = FALSE)
  
  loCRC_matched_ratios[b] <- compute_ratio(g_sub, oral_taxa, crc_taxa_broad, butyrate_producers)
  
  if (b %% 10 == 0) cat("  Completed", b, "/ 100 | ratio =", round(loCRC_matched_ratios[b], 2), "\n")
}

eoCRC_loso <- loso_ratios

cat("\n── Distribution comparison ──\n")
cat("eoCRC LOSO: median =", round(median(eoCRC_loso), 1),
    " IQR =", round(quantile(eoCRC_loso, 0.25), 1), "-",
    round(quantile(eoCRC_loso, 0.75), 1), "\n")
cat("loCRC matched: median =", round(median(loCRC_matched_ratios), 1),
    " IQR =", round(quantile(loCRC_matched_ratios, 0.25), 1), "-",
    round(quantile(loCRC_matched_ratios, 0.75), 1), "\n")

wt <- wilcox.test(eoCRC_loso, loCRC_matched_ratios)
cat("Wilcoxon p =", signif(wt$p.value, 4), "\n")
cat("Cohen's d =", round((mean(eoCRC_loso) - mean(loCRC_matched_ratios)) /
                          sd(loCRC_matched_ratios), 2), "\n")
cat("loCRC exceeding eoCRC median:", sum(loCRC_matched_ratios >= median(eoCRC_loso)),
    "of 100\n")

saveRDS(list(eoCRC_loso = eoCRC_loso, loCRC_matched = loCRC_matched_ratios),
        file.path(NETWORK_DIR, "matched_comparison.rds"))

# ══════════════════════════════════════════════════════════════════
# 6. COMPOSITIONAL HOMOGENEITY
# ══════════════════════════════════════════════════════════════════

cat("\n══ Compositional homogeneity test ══\n\n")

all_crc_idx <- meta_all$group == "CRC"
bc_all_crc <- vegdist(sp[all_crc_idx, ], method = "bray")
age_group <- ifelse(meta_all$age[all_crc_idx] < 50, "eoCRC", "loCRC")

bd <- betadisper(bc_all_crc, age_group)
perm_bd <- permutest(bd, pairwise = TRUE, permutations = 999)

cat("eoCRC distance to centroid:", round(mean(bd$distances[age_group == "eoCRC"]), 4), "\n")
cat("loCRC distance to centroid:", round(mean(bd$distances[age_group == "loCRC"]), 4), "\n")
cat("Betadisper p:", perm_bd$tab$`Pr(>F)`[1], "\n")

saveRDS(list(
  betadisper_eoCRC = mean(bd$distances[age_group == "eoCRC"]),
  betadisper_loCRC = mean(bd$distances[age_group == "loCRC"]),
  betadisper_p = perm_bd$tab$`Pr(>F)`[1]
), file.path(NETWORK_DIR, "compositional_homogeneity.rds"))

# ══════════════════════════════════════════════════════════════════
# 7. HOMOGENEITY-MATCHED TOPOLOGY TEST (30 iterations)
# ══════════════════════════════════════════════════════════════════

cat("\n══ Homogeneity-matched topology test ══\n")
cat("Sampling homogeneous loCRC subsets to test if homogeneity\n")
cat("alone explains the eoCRC topology...\n\n")

bc_lo_full <- as.matrix(vegdist(sp[loCRC_idx_all, ], method = "bray"))

set.seed(42)
n_homo_iter <- 30
homogeneous_ratios <- numeric(n_homo_iter)
random_ratios <- numeric(n_homo_iter)

for (b in 1:n_homo_iter) {
  # Homogeneous subset: pick seed, take 78 most similar
  seed <- sample(1:length(loCRC_idx_all), 1)
  distances_to_seed <- bc_lo_full[seed, ]
  homo_subset <- loCRC_idx_all[order(distances_to_seed)[1:78]]
  
  se_homo <- spiec.easi(sp_filt[homo_subset, ], method = "mb",
                         lambda.min.ratio = 1e-2, nlambda = 20,
                         pulsar.params = list(rep.num = 50, seed = 42))
  adj_homo <- as.matrix(getRefit(se_homo))
  colnames(adj_homo) <- rownames(adj_homo) <- short_names
  g_homo <- graph_from_adjacency_matrix(adj_homo, mode = "undirected", diag = FALSE)
  homogeneous_ratios[b] <- compute_ratio(g_homo, oral_taxa, crc_taxa_broad, butyrate_producers)
  
  # Random subset for comparison
  rand_subset <- sample(loCRC_idx_all, 78)
  se_rand <- spiec.easi(sp_filt[rand_subset, ], method = "mb",
                         lambda.min.ratio = 1e-2, nlambda = 20,
                         pulsar.params = list(rep.num = 50, seed = 42))
  adj_rand <- as.matrix(getRefit(se_rand))
  colnames(adj_rand) <- rownames(adj_rand) <- short_names
  g_rand <- graph_from_adjacency_matrix(adj_rand, mode = "undirected", diag = FALSE)
  random_ratios[b] <- compute_ratio(g_rand, oral_taxa, crc_taxa_broad, butyrate_producers)
  
  if (b %% 5 == 0) cat("  iter", b, ": homogeneous =", round(homogeneous_ratios[b], 1),
                        " random =", round(random_ratios[b], 1), "\n")
}

wt_homo_vs_rand <- wilcox.test(homogeneous_ratios, random_ratios)
wt_homo_vs_eo <- wilcox.test(homogeneous_ratios, eoCRC_loso)

cat("\nHomogeneous loCRC: median =", round(median(homogeneous_ratios), 1), "\n")
cat("Random loCRC: median =", round(median(random_ratios), 1), "\n")
cat("eoCRC LOSO: median =", round(median(eoCRC_loso), 1), "\n")
cat("Homogeneous vs random p =", signif(wt_homo_vs_rand$p.value, 3), "\n")
cat("Homogeneous vs eoCRC p =", signif(wt_homo_vs_eo$p.value, 3), "\n")

saveRDS(list(
  homogeneous_ratios = homogeneous_ratios,
  random_ratios = random_ratios,
  eoCRC_loso = eoCRC_loso,
  p_homo_vs_rand = wt_homo_vs_rand$p.value,
  p_homo_vs_eo = wt_homo_vs_eo$p.value
), file.path(NETWORK_DIR, "homogeneity_topology_test.rds"))

# ══════════════════════════════════════════════════════════════════
# 8. SAVE ALL ROBUSTNESS RESULTS
# ══════════════════════════════════════════════════════════════════

saveRDS(list(
  loso = list(ratios = loso_ratios, studies = loso_names),
  algorithms = algo_results,
  species_sensitivity = sensitivity_results,
  restart_sensitivity = restart_results,
  matched_comparison = list(eoCRC = eoCRC_loso, loCRC = loCRC_matched_ratios),
  homogeneity = list(homogeneous = homogeneous_ratios, random = random_ratios)
), file.path(NETWORK_DIR, "robustness_complete.rds"))

cat("\n══ Step 7 complete. All robustness results saved. ══\n")
cat("Files in", NETWORK_DIR, ":\n")
print(list.files(NETWORK_DIR))
