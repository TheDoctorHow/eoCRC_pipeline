# mechanism/scripts/02_build_network.R
#
# Assemble a unified, Cytoscape-ready PPI network from:
#   - mechanism/data/seed_proteins.tsv         — full seed metadata
#   - mechanism/results/string_nodes_human.tsv — human STRING nodes (seed + shell)
#   - mechanism/results/string_edges_human.tsv — human STRING edges
#   - mechanism/results/string_nodes_p_stomatis.tsv — P. stomatis STRING nodes (may be empty)
#   - mechanism/results/string_edges_p_stomatis.tsv — P. stomatis STRING edges (may be empty)
#   - mechanism/data/bridge_edges.tsv          — literature-curated bacterial→host edges
#   - maaslin2_results_v3/eoCRC_vs_YoungCtrl_species/all_results.tsv — species-level DA results
#
# Outputs:
#   mechanism/results/network_nodes.tsv
#   mechanism/results/network_edges.tsv
#   mechanism/results/network.rds             — igraph object
#
# Run from repo root:
#   Rscript mechanism/scripts/02_build_network.R

# ---------------------------------------------------------------------------
# 0. Dependencies
# ---------------------------------------------------------------------------
for (pkg in c("igraph", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------------------
args      <- commandArgs(trailingOnly = FALSE)
flag_idx  <- grep("^--file=", args)
repo_root <- if (length(flag_idx) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", args[flag_idx])), "../.."))
} else {
  getwd()
}

path_seeds         <- file.path(repo_root, "mechanism", "data", "seed_proteins.tsv")
path_str_nodes_h   <- file.path(repo_root, "mechanism", "results", "string_nodes_human.tsv")
path_str_edges_h   <- file.path(repo_root, "mechanism", "results", "string_edges_human.tsv")
path_str_nodes_ps  <- file.path(repo_root, "mechanism", "results", "string_nodes_p_stomatis.tsv")
path_str_edges_ps  <- file.path(repo_root, "mechanism", "results", "string_edges_p_stomatis.tsv")
path_bridge        <- file.path(repo_root, "mechanism", "data", "bridge_edges.tsv")
path_maaslin2      <- file.path(repo_root, "maaslin2_results_v3",
                                "eoCRC_vs_YoungCtrl_species", "all_results.tsv")
path_str_cache     <- file.path(repo_root, "mechanism", "results",
                                "string_cache_9606", "9606.protein.info.v11.5.txt.gz")
results_dir        <- file.path(repo_root, "mechanism", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 2. Helper: read TSV stripping comment lines
# ---------------------------------------------------------------------------
read_tsv_clean <- function(path, ...) {
  read.table(path, sep = "\t", header = TRUE, comment.char = "#",
             stringsAsFactors = FALSE, quote = "", ...)
}

# ---------------------------------------------------------------------------
# 3. Load sources
# ---------------------------------------------------------------------------

# Seed metadata (human + bacterial; strip comment rows automatically)
seeds <- read_tsv_clean(path_seeds)
cat(sprintf("Seeds loaded: %d rows\n", nrow(seeds)))

# Human STRING nodes / edges
str_nodes_h <- read_tsv_clean(path_str_nodes_h)
str_edges_h <- read_tsv_clean(path_str_edges_h)

# P. stomatis STRING nodes / edges (files exist but may be header-only)
str_nodes_ps <- if (file.exists(path_str_nodes_ps)) {
  df <- read_tsv_clean(path_str_nodes_ps)
  if (nrow(df) == 0L) {
    cat("P. stomatis STRING nodes: file present but empty (taxid 341694 not in STRING v11.5).\n")
  }
  df
} else {
  cat("P. stomatis STRING node file absent — skipping.\n")
  data.frame(STRING_id = character(), node_type = character(), protein = character())
}

str_edges_ps <- if (file.exists(path_str_edges_ps)) {
  df <- read_tsv_clean(path_str_edges_ps)
  if (nrow(df) == 0L) {
    cat("P. stomatis STRING edges: file present but empty.\n")
  }
  df
} else {
  cat("P. stomatis STRING edge file absent — skipping.\n")
  data.frame(from = character(), to = character(), combined_score = numeric())
}

# Bridge edges
bridge <- read_tsv_clean(path_bridge)
cat(sprintf("Bridge edges loaded: %d rows\n", nrow(bridge)))

# MaAsLin2 results — extract bacterial species of interest
maaslin2_raw <- read_tsv_clean(path_maaslin2)
bacterial_species_map <- c(
  "Parvimonas micra"            = "species.Parvimonas.micra",
  "Peptostreptococcus stomatis" = "species.Peptostreptococcus.stomatis"
)
maaslin2 <- maaslin2_raw %>%
  filter(feature %in% bacterial_species_map) %>%
  select(feature, coef, qval)
cat(sprintf("MaAsLin2 bacterial rows found: %d\n", nrow(maaslin2)))

# STRING human protein.info — for resolving shell node gene symbols
str_info_h <- local({
  tmp <- tempfile()
  file.copy(path_str_cache, tmp)
  df  <- read.table(gzfile(tmp), sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, quote = "")
  unlink(tmp)
  # Normalise column names: first col is STRING id, second is preferred_name
  colnames(df)[1:2] <- c("STRING_id", "preferred_name")
  df[, c("STRING_id", "preferred_name")]
})
cat(sprintf("STRING protein.info loaded: %d entries\n", nrow(str_info_h)))

# ---------------------------------------------------------------------------
# 4. Build node table
# ---------------------------------------------------------------------------

# ---- 4a. Human seed nodes ------------------------------------------------
human_seeds <- seeds %>%
  filter(organism == "Homo sapiens") %>%
  transmute(
    node_id         = protein,
    organism        = organism,
    node_type       = "seed_human",
    role            = role,
    uniprot_id      = uniprot_id,
    ncbi_protein_id = as.character(ncbi_protein_id),
    source_pmid     = as.character(source_pmid),
    maaslin2_coef   = NA_real_,
    maaslin2_qval   = NA_real_
  )

# ---- 4b. Human shell nodes -----------------------------------------------
# Join STRING IDs to gene symbols via protein.info cache
shell_str <- str_nodes_h %>%
  filter(node_type == "shell") %>%
  left_join(str_info_h, by = "STRING_id")

unresolved_shell <- sum(is.na(shell_str$preferred_name))
if (unresolved_shell > 0L) {
  warning(sprintf("%d shell nodes could not be resolved to gene symbols.", unresolved_shell))
}

human_shell <- shell_str %>%
  transmute(
    node_id        = ifelse(!is.na(preferred_name), preferred_name, STRING_id),
    organism       = "Homo sapiens",
    node_type      = "shell_human",
    role           = NA_character_,
    uniprot_id     = NA_character_,
    ncbi_protein_id = NA_character_,
    source_pmid    = NA_character_,
    maaslin2_coef  = NA_real_,
    maaslin2_qval  = NA_real_
  )

cat(sprintf("Shell nodes resolved: %d / %d have gene symbols\n",
            sum(!is.na(shell_str$preferred_name)), nrow(shell_str)))

# ---- 4c. Bacterial effector nodes ----------------------------------------
bact_seeds <- seeds %>%
  filter(organism != "Homo sapiens") %>%
  transmute(
    node_id         = protein,
    organism        = organism,
    node_type       = "bacterial_effector",
    role            = role,
    uniprot_id      = uniprot_id,
    ncbi_protein_id = as.character(ncbi_protein_id),
    source_pmid     = as.character(source_pmid),
    maaslin2_coef   = NA_real_,
    maaslin2_qval   = NA_real_
  )

# Attach MaAsLin2 coefficients — keyed by organism name
for (i in seq_len(nrow(bact_seeds))) {
  org   <- bact_seeds$organism[i]
  feat  <- bacterial_species_map[org]
  if (!is.na(feat) && feat %in% maaslin2$feature) {
    row_m <- maaslin2[maaslin2$feature == feat, ]
    bact_seeds$maaslin2_coef[i] <- row_m$coef
    bact_seeds$maaslin2_qval[i] <- row_m$qval
  }
}

# ---- 4d. Combine ----------------------------------------------------------
node_tbl <- bind_rows(human_seeds, human_shell, bact_seeds)
# Remove exact duplicates (should be none, but guard against re-runs)
node_tbl <- node_tbl[!duplicated(node_tbl$node_id), ]

cat(sprintf("\nNode table: %d total nodes\n", nrow(node_tbl)))
cat("  By node_type:\n")
print(table(node_tbl$node_type))

# ---------------------------------------------------------------------------
# 5. Build edge table
# ---------------------------------------------------------------------------

# ---- 5a. STRING-to-gene-symbol lookup table (seeds + shell) --------------
# Maps STRING_id -> node_id (gene symbol) for human network
id_to_symbol <- local({
  # Seeds: from str_nodes_h (protein column = gene symbol)
  seed_map <- str_nodes_h %>%
    filter(node_type == "seed", !is.na(protein), protein != "NA") %>%
    select(STRING_id, node_id = protein)
  # Shell: from shell_str (preferred_name column)
  shell_map <- shell_str %>%
    filter(!is.na(preferred_name)) %>%
    select(STRING_id, node_id = preferred_name)
  bind_rows(seed_map, shell_map)
})

# ---- 5b. Human STRING PPI edges (deduplicated) ---------------------------
# STRING stores undirected edges as A→B AND B→A; keep only canonical direction
# (min(from,to) → max(from,to)) to avoid doubling edge count in Cytoscape
str_edges_h_dedup <- str_edges_h %>%
  mutate(
    canon_from = pmin(from, to),
    canon_to   = pmax(from, to)
  ) %>%
  distinct(canon_from, canon_to, .keep_all = TRUE) %>%
  left_join(id_to_symbol, by = c("canon_from" = "STRING_id")) %>%
  rename(source = node_id) %>%
  left_join(id_to_symbol, by = c("canon_to" = "STRING_id")) %>%
  rename(target = node_id) %>%
  filter(!is.na(source), !is.na(target)) %>%
  transmute(
    source        = source,
    target        = target,
    edge_type     = "STRING_PPI",
    weight        = combined_score,
    pmid          = NA_character_,
    evidence_type = NA_character_
  )

cat(sprintf("\nHuman STRING edges: %d raw -> %d after deduplication\n",
            nrow(str_edges_h), nrow(str_edges_h_dedup)))

# ---- 5c. P. stomatis STRING edges (if any) --------------------------------
pstom_edges <- if (nrow(str_edges_ps) > 0L) {
  # Build a symbol map for p. stomatis nodes (FBA has no STRING entry, so
  # this block only executes if STRING ever indexes taxid 341694 in future)
  ps_sym <- if (nrow(str_nodes_ps) > 0L) {
    str_nodes_ps %>%
      filter(!is.na(protein), protein != "NA") %>%
      select(STRING_id, node_id = protein)
  } else {
    data.frame(STRING_id = character(), node_id = character())
  }
  all_sym <- bind_rows(id_to_symbol, ps_sym)

  str_edges_ps %>%
    mutate(canon_from = pmin(from, to), canon_to = pmax(from, to)) %>%
    distinct(canon_from, canon_to, .keep_all = TRUE) %>%
    left_join(all_sym, by = c("canon_from" = "STRING_id")) %>%
    rename(source = node_id) %>%
    left_join(all_sym, by = c("canon_to" = "STRING_id")) %>%
    rename(target = node_id) %>%
    filter(!is.na(source), !is.na(target)) %>%
    transmute(source, target,
              edge_type     = "STRING_PPI",
              weight        = combined_score,
              pmid          = NA_character_,
              evidence_type = NA_character_)
} else {
  cat("No P. stomatis STRING edges to add.\n")
  data.frame(source = character(), target = character(),
             edge_type = character(), weight = numeric(),
             pmid = character(), evidence_type = character())
}

# ---- 5d. Literature-curated bridge edges ---------------------------------
bridge_edges <- bridge %>%
  transmute(
    source        = source_protein,
    target        = target_protein,
    edge_type     = "literature_curated",
    weight        = NA_real_,
    pmid          = as.character(source_pmid),
    evidence_type = evidence_type
  )

cat(sprintf("Bridge edges: %d\n", nrow(bridge_edges)))

# ---- 5e. Combine all edges -----------------------------------------------
edge_tbl <- bind_rows(str_edges_h_dedup, pstom_edges, bridge_edges)

cat(sprintf("\nEdge table: %d total edges\n", nrow(edge_tbl)))
cat("  By edge_type:\n")
print(table(edge_tbl$edge_type))

# ---------------------------------------------------------------------------
# 6. Isolated node check
# ---------------------------------------------------------------------------
connected_nodes <- unique(c(edge_tbl$source, edge_tbl$target))
isolated <- setdiff(node_tbl$node_id, connected_nodes)
if (length(isolated) > 0L) {
  cat(sprintf("\nWARNING: %d isolated node(s) with no edges:\n  %s\n",
              length(isolated), paste(isolated, collapse = ", ")))
} else {
  cat("\nAll nodes have at least one edge.\n")
}

# ---------------------------------------------------------------------------
# 7. Build igraph object
# ---------------------------------------------------------------------------
g <- graph_from_data_frame(
  d        = edge_tbl[, c("source", "target")],
  directed = FALSE,
  vertices = node_tbl$node_id
)

# Attach node attributes
for (col in setdiff(colnames(node_tbl), "node_id")) {
  g <- set_vertex_attr(g, col, value = node_tbl[[col]])
}

# Attach edge attributes
for (col in setdiff(colnames(edge_tbl), c("source", "target"))) {
  g <- set_edge_attr(g, col, value = edge_tbl[[col]])
}

cat(sprintf("\nigraph: %d vertices, %d edges\n", vcount(g), ecount(g)))

# ---------------------------------------------------------------------------
# 8. Save outputs
# ---------------------------------------------------------------------------
node_out <- file.path(results_dir, "network_nodes.tsv")
edge_out  <- file.path(results_dir, "network_edges.tsv")
rds_out   <- file.path(results_dir, "network.rds")

write.table(node_tbl, node_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edge_tbl, edge_out, sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(g, rds_out)

cat(sprintf("Saved:\n  %s\n  %s\n  %s\n", node_out, edge_out, rds_out))

cat("\n==================== SUMMARY ====================\n")
cat(sprintf("Total nodes : %d\n", nrow(node_tbl)))
cat(sprintf("  seed_human        : %d\n", sum(node_tbl$node_type == "seed_human")))
cat(sprintf("  shell_human       : %d\n", sum(node_tbl$node_type == "shell_human")))
cat(sprintf("  bacterial_effector: %d\n", sum(node_tbl$node_type == "bacterial_effector")))
cat(sprintf("Total edges : %d\n", nrow(edge_tbl)))
cat(sprintf("  STRING_PPI         : %d\n", sum(edge_tbl$edge_type == "STRING_PPI")))
cat(sprintf("  literature_curated : %d\n", sum(edge_tbl$edge_type == "literature_curated")))
if (length(isolated) > 0L) {
  cat(sprintf("Isolated nodes: %s\n", paste(isolated, collapse = ", ")))
} else {
  cat("Isolated nodes: none\n")
}
cat("=================================================\n")
