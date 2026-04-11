# Mechanism Module — Bacterial-Host PPI Network

## Goal

Build a protein–protein interaction (PPI) network connecting bacterial effectors from
four eoCRC-enriched anaerobes to human host signalling/epigenetic/tumour-suppressor
proteins, in order to propose mechanistic hypotheses for how these microbes may promote
early-onset colorectal carcinogenesis.

### Target bacterial species
| Abbreviation | Full name |
|---|---|
| P. micra | *Parvimonas micra* |
| P. stomatis | *Parvimonas stomatis* |
| G. morbillorum | *Gemella morbillorum* |
| D. pneumosintes | *Dialister pneumosintes* |

---

## Planned workflow

1. **Seed list curation** (`data/seed_proteins.tsv`)
   Host proteins from oncogenic, inflammatory, epigenetic, and DNA-repair pathways +
   known/predicted bacterial effectors. UniProt IDs and PMIDs verified manually.

2. **PPI retrieval** (`scripts/01_fetch_ppi.R` or `.py`)
   Query STRING (score ≥ 700), BioGRID, and IntAct for all pairwise interactions
   among seed proteins. Also pull bacterial–host cross-kingdom interactions where
   annotated (HPIDB, PHI-base).

3. **Network construction** (`scripts/02_build_network.R`)
   Combine retrieved edges, filter by evidence tier, build igraph object.
   Node attributes: organism, role, log2FC from MaAsLin2 v3 (where available).

4. **Network analysis** (`scripts/03_analyse_network.R`)
   Hub scoring, shortest-path distances from bacterial effectors to tumour suppressors,
   community detection, enrichment of CRC GWAS genes.

5. **Visualisation** (`scripts/04_plot_network.R`)
   Publication-quality Cytoscape-style plot saved to `figures/`.

---

## Directory layout

```
mechanism/
├── README.md           — this file
├── data/
│   └── seed_proteins.tsv   — curated seed list (starting point for PPI queries)
├── scripts/            — analysis scripts (numbered by step)
├── results/            — intermediate tables, edge lists, community tables
└── figures/            — network plots and summary figures
```

---

## Key references

- Flemer et al. 2018 — *P. micra* enrichment in CRC (gut microbiome profiling)
- Tsoi et al. / Koliarakis et al. — bacterial mechanisms in CRC
- STRING database — Szklarczyk et al. 2023
- PHI-base / HPIDB — host–pathogen interaction databases
