#!/usr/bin/env Rscript
# Figure 5 — Key species boxplots (v3 q-values; Faecalibacterium prausnitzii replaces Intestinimonas)

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
  library(dplyr);   library(tidyr)
  library(compositions)
})

if (basename(getwd()) == "scripts") setwd("..")
dir.create("figures", showWarnings = FALSE)

if (!file.exists("all_samples_metadata.rds") ||
    !file.exists("species_filtered.rds"))
  stop("Run scripts/step1_extract_filter.R first to generate input data.")
if (!file.exists("maaslin2_results_v3/summaries/all_associations.csv"))
  stop("Run scripts/step2c_maaslin2_v3.R first to generate MaAsLin2 v3 results.")

# ── Palette & theme ────────────────────────────────────────────────────────────
COL <- list(
  eoCRC      = "#E74C3C", young_ctrl = "#27AE60",
  loCRC      = "#2980B9", older_ctrl = "#8E44AD"
)

theme_pres <- function(base = 12)
  theme_bw(base_size = base) +
  theme(plot.title    = element_text(face = "bold", size = base + 1),
        plot.subtitle = element_text(size = base - 1, color = "grey40"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"))

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading data...\n")
meta    <- readRDS("all_samples_metadata.rds") %>%
  mutate(grp = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))
sp_mat  <- readRDS("species_filtered.rds")

# Load v3 MaAsLin2 results (canonical CLR, fixed)
maas_v3 <- read.csv("maaslin2_results_v3/summaries/all_associations.csv",
                    stringsAsFactors = FALSE)

pseudo <- 1e-6

# CLR-transform full species matrix (for display only — data viz, not ML)
sp_clr_all <- as.matrix(clr(sp_mat + pseudo))
# Column names are already in "species:Name" format — preserve them
colnames(sp_clr_all) <- colnames(sp_mat)

# ── Key species (Faecalibacterium replaces Intestinimonas) ────────────────────
key_sp <- c("species:Parvimonas micra",
            "species:Gemella morbillorum",
            "species:Dialister pneumosintes",
            "species:Faecalibacterium prausnitzii")

# ── Q-value lookup from v3 results ────────────────────────────────────────────
get_q_v3 <- function(sp_raw) {
  k  <- gsub("^species:", "", sp_raw)
  mn <- paste0("species.", gsub(" ", "\\.", k))
  r  <- maas_v3 %>%
    filter(comparison == "eoCRC_vs_YoungCtrl", feature_set == "species",
           feature == mn) %>%
    pull(qval)
  if (length(r) == 0) NA_real_ else r[1]
}
q_annot <- sapply(key_sp, get_q_v3)
names(q_annot) <- key_sp
cat("Q-values (v3):\n"); print(q_annot)

# ── Q-label formatter: scientific notation for q < 0.001, else 3 d.p. ─────────
make_qlab <- function(q) {
  if (is.na(q)) return("ns")
  if (q < 0.001) {
    # Format as q=X.Xe-XX
    sci <- formatC(q, format = "e", digits = 1)
    return(paste0("q=", sci))
  }
  if (q < 0.05)  return(sprintf("q=%.3f", q))
  sprintf("q=%.3f", q)
}

# ── Build boxplot data ─────────────────────────────────────────────────────────
# Check which key species are present in the CLR matrix
avail_sp <- intersect(key_sp, colnames(sp_clr_all))
missing_sp <- setdiff(key_sp, avail_sp)
if (length(missing_sp) > 0) {
  cat("WARNING: these species not found in matrix:\n"); print(missing_sp)
}

box_data <- sp_clr_all[, avail_sp, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(.), grp = meta$grp[match(rownames(.), rownames(sp_clr_all))]) %>%
  filter(!is.na(grp)) %>%
  pivot_longer(-c(sample_id, grp), names_to = "feature", values_to = "clr") %>%
  mutate(
    grp       = factor(grp, levels = c("eoCRC", "young_ctrl", "loCRC", "older_ctrl")),
    grp_label = recode(grp,
                       eoCRC      = "eoCRC",
                       young_ctrl = "Young\nctrl",
                       loCRC      = "loCRC",
                       older_ctrl = "Older\nctrl"),
    grp_label = factor(grp_label,
                       levels = c("eoCRC", "Young\nctrl", "loCRC", "Older\nctrl"))
  )

# Bracket y positions: leave enough room above the data for bar + text
ymax_sp <- box_data %>%
  group_by(feature) %>%
  summarise(ymax = max(clr, na.rm = TRUE) + 0.5, .groups = "drop") %>%
  { setNames(.$ymax, .$feature) }

# ── Per-species titles: italicised + depletion note for Faecalibacterium ───────
sp_titles <- c(
  "species:Parvimonas micra"           = "Parvimonas micra",
  "species:Gemella morbillorum"        = "Gemella morbillorum",
  "species:Dialister pneumosintes"     = "Dialister pneumosintes",
  "species:Faecalibacterium prausnitzii" = "Faecalibacterium prausnitzii"
)

# ── Build panels ──────────────────────────────────────────────────────────────
panel_plots <- lapply(avail_sp, function(sp) {
  d     <- filter(box_data, feature == sp)
  q_txt <- make_qlab(q_annot[sp])
  ytop  <- ymax_sp[sp]
  title <- sp_titles[sp]

  ggplot(d, aes(x = grp_label, y = clr, fill = grp)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, color = "grey30") +
    geom_jitter(aes(color = grp), width = 0.18, alpha = 0.5, size = 0.9) +
    scale_fill_manual(
      values = c("eoCRC"      = COL$eoCRC,
                 "Young\nctrl" = COL$young_ctrl,
                 "loCRC"       = COL$loCRC,
                 "Older\nctrl" = COL$older_ctrl),
      guide = "none"
    ) +
    scale_color_manual(
      values = c("eoCRC"      = COL$eoCRC,
                 "Young\nctrl" = COL$young_ctrl,
                 "loCRC"       = COL$loCRC,
                 "Older\nctrl" = COL$older_ctrl),
      guide = "none"
    ) +
    # Significance bracket: eoCRC (x=1) vs Young ctrl (x=2)
    # horizontal bar at ytop; text sits 0.55 above the bar with clear gap
    annotate("segment", x = 1,   xend = 2,   y = ytop,        yend = ytop,        linewidth = 0.5) +
    annotate("segment", x = 1,   xend = 1,   y = ytop - 0.2,  yend = ytop,        linewidth = 0.5) +
    annotate("segment", x = 2,   xend = 2,   y = ytop - 0.2,  yend = ytop,        linewidth = 0.5) +
    annotate("text",    x = 1.5,             y = ytop + 0.55, label = q_txt, size = 3) +
    labs(title    = title,
         subtitle = if (sp == "species:Faecalibacterium prausnitzii")
                      "Butyrate producer depleted in eoCRC (q=0.082)" else NULL,
         x = NULL, y = "CLR-transformed abundance") +
    theme_pres(10) +
    theme(plot.title = element_text(face = "bold.italic", size = 9))
})

# ── Compose & save ─────────────────────────────────────────────────────────────
fig5 <- wrap_plots(panel_plots, ncol = 2) +
  plot_annotation(
    title    = "Key species abundance across groups \u2014 CLR-transformed",
    subtitle = paste0("Boxes: IQR | Points: individual samples | ",
                      "Brackets: MaAsLin2 v3 q-value (eoCRC vs Young ctrl)"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("figures/fig5_species_boxplots.png", fig5, width = 11, height = 8, dpi = 300,
       device = "png")
ggsave("figures/fig5_species_boxplots.pdf", fig5, width = 11, height = 8,
       device = cairo_pdf)
cat("Saved: figures/fig5_species_boxplots.png + .pdf\n")
cat("Done.\n")
