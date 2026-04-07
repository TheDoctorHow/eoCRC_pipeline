#!/usr/bin/env Rscript
# Figure 7 — Biological narrative: clean two-column text comparison

suppressPackageStartupMessages({ library(ggplot2) })

if (basename(getwd()) == "scripts") setwd("..")
dir.create("figures", showWarnings = FALSE)

# ── Size constants (ggplot2 annotate: 1 pt ≈ 1/2.845 size units) ─────────────
SZ_TITLE   <- 7.2   # ≈20 pt
SZ_HEADER  <- 6.4   # ≈18 pt
SZ_ITEM    <- 5.0   # ≈14 pt
SZ_SUB     <- 3.9   # ≈11 pt  (parenthetical sub-labels)
SZ_SUMMARY <- 4.2   # ≈12 pt
SZ_CAPTION <- 3.2   # ≈ 9 pt

# ── Colour scheme ─────────────────────────────────────────────────────────────
GRN_HD  <- "#1E8449"   # header text green
GRN_BG  <- "#D5F5E3"   # header rect green
GRN_CHK <- "#1E8449"   # checkmark colour
RED_HD  <- "#922B21"   # header text red
RED_BG  <- "#FADBD8"   # header rect red
RED_X   <- "#C0392B"   # ✗ colour
ARROW_C <- "#2C3E50"   # centre arrow
BOX_BG  <- "#F2F3F4"   # summary box background

# ── Coordinate space: x 0-10, y 0-10 ─────────────────────────────────────────
# Left col x: 0.15 – 4.50  |  Centre: 4.50 – 5.50  |  Right col: 5.50 – 9.85
L_MARK  <- 0.35   # x of ✓ / ✗ marks
L_TEXT  <- 0.90   # x of left item text (hjust=0)
R_MARK  <- 5.72   # x of ✗ marks (right column)
R_TEXT  <- 6.27   # x of right item text (hjust=0)
CTR     <- 5.00   # centre of arrow
L_HEAD  <- 2.30   # centre of left header
R_HEAD  <- 7.68   # centre of right header

# Item y positions — right column has multi-line items so needs more spacing
# 6 items each side; usable y: 8.35 (below header) down to 2.30 (above summary)
# right col item centres (allow extra gap for multi-line items)
R_Y <- c(7.90, 6.90, 5.95, 5.20, 4.50, 3.55)
# left col items align to same rows
L_Y <- R_Y

# Arrow placed at vertical midpoint of item range
ARR_Y <- mean(R_Y)   # ≈5.92

fig7 <- ggplot() +

  # ── Title & subtitle ─────────────────────────────────────────────────────────
  annotate("text", x=CTR, y=9.78,
           label="eoCRC Microbiome Signature \u2014 Proposed Mechanism",
           size=SZ_TITLE, fontface="bold", color="grey10", hjust=0.5) +
  annotate("text", x=CTR, y=9.35,
           label="Oral-gut axis hypothesis",
           size=SZ_CAPTION+0.8, color="grey45", hjust=0.5, fontface="italic") +

  # ── Vertical separator lines — drawn FIRST so header rects cover them ─────────
  # Stop at y=8.45 (below header bottom at 8.55) to avoid cutting header visually
  annotate("segment", x=4.50, xend=4.50, y=2.25, yend=8.45,
           color="grey80", linewidth=0.5, linetype="dashed") +
  annotate("segment", x=5.50, xend=5.50, y=2.25, yend=8.45,
           color="grey80", linewidth=0.5, linetype="dashed") +

  # ── Column header backgrounds ─────────────────────────────────────────────────
  # Green extends from x=0 (full left edge) to x=4.50 (at separator)
  annotate("rect", xmin=0.00, xmax=4.50, ymin=8.55, ymax=9.05,
           fill=GRN_BG, color=NA) +
  annotate("rect", xmin=5.50, xmax=10.00, ymin=8.55, ymax=9.05,
           fill=RED_BG, color=NA) +

  # ── Column headers ────────────────────────────────────────────────────────────
  annotate("text", x=L_HEAD, y=8.80,
           label="\u2713 HEALTHY YOUNG GUT",
           size=SZ_HEADER, fontface="bold", color=GRN_HD, hjust=0.5) +
  annotate("text", x=R_HEAD, y=8.80,
           label="\u2717 eoCRC GUT",
           size=SZ_HEADER, fontface="bold", color=RED_HD, hjust=0.5) +

  # ── Thin divider lines under headers (full column width) ─────────────────────
  annotate("segment", x=0.00, xend=4.50, y=8.55, yend=8.55,
           color=GRN_HD, linewidth=0.4) +
  annotate("segment", x=5.50, xend=10.00, y=8.55, yend=8.55,
           color=RED_HD, linewidth=0.4) +

  # ════════════════════════════════════════════════════════════════════════════
  # LEFT COLUMN — Healthy items
  # ════════════════════════════════════════════════════════════════════════════

  # Item 1 — species name italic, descriptor normal
  annotate("text", x=L_MARK, y=L_Y[1], label="\u2713",
           size=SZ_ITEM, color=GRN_CHK, hjust=0.5, fontface="bold") +
  annotate("text", x=L_TEXT, y=L_Y[1],
           label="F. prausnitzii abundant",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 2
  annotate("text", x=L_MARK, y=L_Y[2], label="\u2713",
           size=SZ_ITEM, color=GRN_CHK, hjust=0.5, fontface="bold") +
  annotate("text", x=L_TEXT, y=L_Y[2],
           label="Butyrate production high",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 3
  annotate("text", x=L_MARK, y=L_Y[3], label="\u2713",
           size=SZ_ITEM, color=GRN_CHK, hjust=0.5, fontface="bold") +
  annotate("text", x=L_TEXT, y=L_Y[3],
           label="Tight junctions intact",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 4
  annotate("text", x=L_MARK, y=L_Y[4], label="\u2713",
           size=SZ_ITEM, color=GRN_CHK, hjust=0.5, fontface="bold") +
  annotate("text", x=L_TEXT, y=L_Y[4],
           label="Colonization resistance maintained",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 5
  annotate("text", x=L_MARK, y=L_Y[5], label="\u2713",
           size=SZ_ITEM, color=GRN_CHK, hjust=0.5, fontface="bold") +
  annotate("text", x=L_TEXT, y=L_Y[5],
           label="Tryptophan metabolism active",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 6
  annotate("text", x=L_MARK, y=L_Y[6], label="\u2713",
           size=SZ_ITEM, color=GRN_CHK, hjust=0.5, fontface="bold") +
  annotate("text", x=L_TEXT, y=L_Y[6],
           label="Low oral pathobiont abundance",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # ════════════════════════════════════════════════════════════════════════════
  # RIGHT COLUMN — eoCRC items
  # ════════════════════════════════════════════════════════════════════════════

  # Item 1 (multi-line: main + parenthetical species)
  annotate("text", x=R_MARK, y=R_Y[1], label="\u2717",
           size=SZ_ITEM, color=RED_X, hjust=0.5, fontface="bold") +
  annotate("text", x=R_TEXT, y=R_Y[1]+0.22,
           label="Oral pathogens colonize colon",
           size=SZ_ITEM, color="grey15", hjust=0) +
  annotate("text", x=R_TEXT, y=R_Y[1]-0.30,
           label="(Parvimonas micra,  Gemella morbillorum,\nPeptostreptococcus stomatis)",
           size=SZ_SUB, color="grey35", hjust=0, fontface="italic") +

  # Item 2 (multi-line: main + species)
  annotate("text", x=R_MARK, y=R_Y[2], label="\u2717",
           size=SZ_ITEM, color=RED_X, hjust=0.5, fontface="bold") +
  annotate("text", x=R_TEXT, y=R_Y[2]+0.20,
           label="Butyrate producers depleted",
           size=SZ_ITEM, color="grey15", hjust=0) +
  annotate("text", x=R_TEXT, y=R_Y[2]-0.25,
           label="(F. prausnitzii\u2193,  Intestinimonas\u2193)",
           size=SZ_SUB, color="grey35", hjust=0, fontface="italic") +

  # Item 3
  annotate("text", x=R_MARK, y=R_Y[3], label="\u2717",
           size=SZ_ITEM, color=RED_X, hjust=0.5, fontface="bold") +
  annotate("text", x=R_TEXT, y=R_Y[3],
           label="Leaky gut barrier",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 4
  annotate("text", x=R_MARK, y=R_Y[4], label="\u2717",
           size=SZ_ITEM, color=RED_X, hjust=0.5, fontface="bold") +
  annotate("text", x=R_TEXT, y=R_Y[4],
           label="Colonization resistance failure",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 5
  annotate("text", x=R_MARK, y=R_Y[5], label="\u2717",
           size=SZ_ITEM, color=RED_X, hjust=0.5, fontface="bold") +
  annotate("text", x=R_TEXT, y=R_Y[5],
           label="Tryptophan biosynthesis depleted",
           size=SZ_ITEM, color="grey15", hjust=0) +

  # Item 6 (multi-line)
  annotate("text", x=R_MARK, y=R_Y[6], label="\u2717",
           size=SZ_ITEM, color=RED_X, hjust=0.5, fontface="bold") +
  annotate("text", x=R_TEXT, y=R_Y[6]+0.20,
           label="Inflammatory environment",
           size=SZ_ITEM, color="grey15", hjust=0) +
  annotate("text", x=R_TEXT, y=R_Y[6]-0.22,
           label="(LPS + FadA invasion)",
           size=SZ_SUB, color="grey40", hjust=0, fontface="italic") +

  # ── Centre dysbiosis arrow ────────────────────────────────────────────────────
  annotate("segment",
           x=4.60, xend=5.40, y=ARR_Y, yend=ARR_Y,
           arrow=arrow(length=unit(0.35,"cm"), type="closed"),
           linewidth=2.8, color=ARROW_C) +
  annotate("text", x=CTR, y=ARR_Y+0.55,
           label="Dysbiosis", size=4.4, fontface="bold", color=ARROW_C, hjust=0.5) +

  # ── Summary box ───────────────────────────────────────────────────────────────
  annotate("rect", xmin=0.15, xmax=9.85, ymin=1.40, ymax=2.15,
           fill=BOX_BG, color=NA) +
  annotate("text", x=CTR, y=1.78,
           label=paste0("eoCRC signature: enriched oral pathobionts + depleted butyrate producers\n",
                        "\u2192 weakened mucosal protection \u2192 malignant transformation"),
           size=SZ_SUMMARY, fontface="bold", color="grey15", hjust=0.5, lineheight=1.3) +

  # ── Caption ───────────────────────────────────────────────────────────────────
  annotate("text", x=CTR, y=1.15,
           label=paste0("Proposed mechanism \u2014 associative evidence only;",
                        " causal direction requires prospective validation"),
           size=SZ_CAPTION, color="grey55", hjust=0.5, fontface="italic") +

  # ── Coordinate / theme ───────────────────────────────────────────────────────
  coord_cartesian(xlim=c(0,10), ylim=c(0.6, 10), expand=FALSE) +
  labs(title=NULL) +
  theme_void() +
  theme(plot.background = element_rect(fill="white", color=NA),
        plot.margin     = margin(4, 6, 4, 6))

ggsave("figures/fig7_biological_narrative.png", fig7,
       width=12, height=7, dpi=300, bg="white")
ggsave("figures/fig7_biological_narrative.pdf", fig7,
       width=12, height=7, device=cairo_pdf)
cat("Saved: figures/fig7_biological_narrative.png + .pdf\n")
cat("Done.\n")
