# eoCRC Microbiome Network Analysis — Technical Report

**Date:** April 4, 2026  
**Dataset:** 9-cohort meta-analysis, n=1,282 (78 eoCRC, 121 young controls, 562 loCRC, 521 older controls)  
**Data source:** curatedMetagenomicData v3, MetaPhlAn4 species profiles  
**Network method:** SpiecEasi (MB method, StARS selection, rep.num=50)

---

## 1. Dual-Module Architecture of CRC Dysbiosis

### Finding
CRC is characterized by two anticorrelated ecological modules that occupy completely separate network communities across all four groups:

**Oral-CRC pathobiont module:** P. micra, G. morbillorum, P. stomatis, B. fragilis, F. nucleatum, H. hathewayi  
**Butyrate producer module:** F. prausnitzii, R. intestinalis, R. inulinivorans, A. hadrus, C. eutactus

### Evidence
- The two modules are anticorrelated in every group: eoCRC rho=-0.377 (p=7.4e-4), young controls rho=-0.268 (p=0.003), loCRC rho=-0.294 (p=1.4e-12), older controls rho=-0.371 (p<1e-16).
- Zero shared community membership between modules in eoCRC and young controls (Louvain detection). One shared community in loCRC (H. hathewayi and F. prausnitzii both in community 1), but this reflects H. hathewayi being a gut-resident species, not true module overlap.
- Zero direct network edges between any butyrate producer and any oral-CRC species in any group.
- The dual-module ratio (oral-CRC score minus butyrate score) predicts CRC status after study adjustment (logistic regression p=7.4e-25 for all CRC vs controls; p=0.003 for eoCRC vs young controls).
- Dual ratio AUC=0.673 for all CRC vs controls, outperforming any single species (F. prausnitzii AUC=0.604, P. micra AUC=0.644).

### Assessment: ROBUST FINDING
This is consistent across all groups, survives study adjustment, and has clear biological interpretation. It is not eoCRC-specific — it characterizes CRC generally.

### Limitation
The dual ratio does not distinguish eoCRC from loCRC (p=0.651). The module species are present at the same abundance and co-occurrence rates in both CRC groups.

---

## 2. Network Diffusion: Directional Ecological Influence Flow

### Finding
When a random walk with restart is seeded at oral pathobiont nodes (P. micra, G. morbillorum, P. stomatis, D. pneumosintes), the diffusion reaches CRC-driver species vs butyrate producers at dramatically different ratios depending on the group:

| Group | CRC:BP Diffusion Ratio | Permutation p |
|-------|----------------------|---------------|
| eoCRC | 40.4 | <0.001 |
| Young controls | 2.9 | 0.31 |
| loCRC (full, n=562) | 2.3 | — |
| Older controls | [computed] | — |
| Permutation mean | 2.5 | — |
| Permutation 99th %ile | 20.4 | — |

In eoCRC, the top species reached by diffusion from oral pathobionts are CRC drivers: E. clostridioformis (influence=0.089), C. symbiosum (0.022), C. innocuum (0.021), E. tayi (0.017), H. hathewayi (0.016), B. fragilis (0.015). In controls, the top species reached are oral commensals: Eubacterium sulci, Veillonella, Streptococcus, Gemella haemolysans.

### Robustness checks

**Leave-one-study-out (9 iterations):** Ratio ranged from 8.8 (drop YachidaS) to 25.7 (drop ZellerG). All iterations remained well above control baseline of 2-3. No single study drives the result.

**Community detection algorithm sensitivity:** Louvain (mean 8 co-clusterings), walktrap (10), infomap (7), label propagation (7.3) — all show eoCRC co-clustering at 100% of runs. Young controls show 0.1-4 depending on algorithm.

**Restart probability sensitivity:** eoCRC consistently exceeds controls and loCRC at all restart values (0.05-0.50). The absolute ratio changes but relative ordering is invariant.

**Species definition sensitivity:** Tested strict (3 CRC species, 2 BP), broad (7 CRC, 5 BP), and very broad definitions. eoCRC-to-control fold difference ranged 8.7x-14.5x regardless of definition.

### CRITICAL TEST: Sample-size-matched loCRC comparison (100 iterations)

To address the sample size asymmetry (eoCRC n=78 vs loCRC n=562), we subsampled loCRC to n=78 across 100 independent iterations and compared the resulting diffusion ratio distribution against the eoCRC leave-one-study-out distribution (9 values).

| Metric | eoCRC (LOSO, 9 values) | loCRC matched n=78 (100 values) |
|--------|----------------------|-------------------------------|
| Mean | 17.2 | 10.3 |
| Median | 16.0 | 5.4 |
| IQR | 14.4 - 19.9 | 2.2 - 10.7 |
| Range | 8.8 - 25.7 | 0.2 - 88.7 |
| CV | 0.33 | 1.50 |

**Wilcoxon rank-sum test: p = 0.001**  
**Cohen's d: 0.46 (moderate effect)**

Key overlap statistics: 83% of loCRC subsamples fall below the eoCRC median (16.0). Only 17 of 100 loCRC subsamples exceed the eoCRC median. Only 29 of 100 exceed even the eoCRC minimum (8.8).

The eoCRC distribution is remarkably tight (CV=0.33) — every leave-one-study-out iteration produces a ratio between 8.8 and 25.7, regardless of which study is dropped. The loCRC distribution is highly variable (CV=1.50) — at n=78, which 78 loCRC patients are sampled dramatically changes the network topology. This variability itself is informative: it suggests loCRC does not have a single characteristic network topology, while eoCRC does.

### Assessment: NOVEL, STATISTICALLY SIGNIFICANT, MODERATE EFFECT

The eoCRC-loCRC topological distinction is statistically real (p=0.001) but moderate in effect size (d=0.46). The eoCRC network consistently produces high oral-to-CRC diffusion ratios across all study compositions, while loCRC network topology is unstable and usually lower but occasionally matches eoCRC by chance. The most defensible interpretation is:

"eoCRC harbors a more stereotyped microbial ecological architecture than loCRC, characterized by consistently elevated oral-to-CRC influence flow (LOSO median 16.0 vs matched loCRC median 5.4, p=0.001). The high topological variability in loCRC (CV=1.50 vs eoCRC CV=0.33) suggests that late-onset CRC encompasses diverse microbial ecologies, while early-onset CRC converges on a more uniform network structure in which oral pathobionts are reliably wired to CRC-driver species."

---

## 3. Within-Study Reproducibility of Key Co-occurrences

### Finding
Pairwise Spearman correlations among eoCRC samples were tested within individual cohorts:

**Gemella morbillorum ↔ P. stomatis:** 5/5 studies positive. Significant within WirbelJ (rho=0.885, p=0.008) and YachidaS (rho=0.677, p=2.1e-5). MOST REPRODUCIBLE.

**P. micra ↔ F. nucleatum:** All testable studies positive. Significant within WirbelJ (rho=0.832, p=0.02) and YachidaS (rho=0.526, p=0.002). STRONG.

**P. stomatis ↔ B. fragilis:** 4/5 positive, but near zero in YachidaS (rho=0.015). Pooled rho=0.324, p=0.004. Partial correlation controlling for study: rho=0.323, p=0.004. POOLED SIGNAL ROBUST, WITHIN-STUDY INCONSISTENT.

**P. micra ↔ B. fragilis:** Mixed directions across studies. NOT REPRODUCIBLE at pairwise level.

**P. stomatis ↔ F. nucleatum:** All testable studies positive. Pooled rho=0.264, p=0.02. DIRECTIONALLY CONSISTENT.

### Assessment: PARTIALLY REPRODUCIBLE
The oral-oral partnerships (Gemella-P. stomatis, P. micra-F. nucleatum) replicate robustly within individual studies. The oral-to-CRC-driver connections (P. stomatis-B. fragilis) are significant in pooled analysis and survive study adjustment but do not replicate within the largest individual cohort (YachidaS).

---

## 4. Community Detection: Oral-CRC Co-clustering

### Finding
Louvain community detection revealed that oral pathobionts co-cluster with CRC-associated species specifically in eoCRC:

| Group | Co-clusterings | Permutation p | Stable across 100 runs |
|-------|---------------|---------------|----------------------|
| eoCRC | 7 | 0.01 | 100% (mean 7.6, range 6-12) |
| Young controls | 0 | — | 32% of runs show >0 |
| loCRC | 0 | — | 0% of runs |
| Older controls | 0 | — | 54% of runs show >0 |

Specific co-clustering in eoCRC: P. stomatis with B. fragilis + H. hathewayi + C. symbiosum (community 14); P. micra + D. pneumosintes with E. clostridioformis (community 23); Gemella with C. innocuum + E. tayi (community 5).

### Robustness
- Leave-one-study-out: 100% positive in all 9 drops (range 4.9-11.8 mean co-clusterings)
- Algorithm robustness: Louvain, walktrap, infomap, label propagation all show eoCRC signal
- Sample-size-matched loCRC: mean 2.6 co-clusterings (50% positive), vs eoCRC 7-8

### CRITICAL CAVEAT
Per-sample scores (co-presence, co-abundance, integration ratio) do not distinguish eoCRC from loCRC (all p>0.4). The community detection finding is a group-level network property that does not translate into measurable per-sample differences. This could reflect genuine higher-order ecological structure visible only at the network level, or it could reflect estimation artifacts from comparing networks built at different sample sizes.

### Assessment: SIGNIFICANT BUT AMBIGUOUS
The result is statistically significant, robust across algorithms and studies, but its biological interpretation is uncertain because it doesn't manifest at the per-sample level.

---

## 5. Patient Stratification: Integrated vs Non-Integrated

### Finding
Defining "integrated" patients as those with ≥2 oral pathobionts AND ≥2 CRC drivers present:

| Group | % Integrated | 
|-------|-------------|
| eoCRC | 25.6% (20/78) |
| loCRC | 21.9% (123/562) |
| Fisher's p | 0.47 |

The proportion is indistinguishable between eoCRC and loCRC.

However, within eoCRC, integrated patients show:
- Lower butyrate producer score (1.21 vs 1.91, Wilcoxon p=0.034)
- More extreme dual-module ratio (-0.44 vs -2.44, p=2.4e-5)
- Elevated functional pathway activity for every CRC-driver species: B. fragilis 2.9x (p=1.9e-4), H. hathewayi 4.3x (p=0.005), P. micra 13.4x (p=4.3e-8), F. nucleatum 5.2x (p=1.3e-6)
- Trending higher species richness (97.3 vs 80.7, p=0.065)

### CRITICAL RESULT
The same functional pattern appears in loCRC integrated patients: B. fragilis 2.2x (p=4.9e-19), P. micra 15x (p=1.2e-38), F. nucleatum 13x (p=1.6e-25). The functional elevation in integrated patients is NOT eoCRC-specific.

### Assessment: REAL BIOLOGY, NOT eoCRC-SPECIFIC
The integrated patient subgroup is biologically meaningful — they have the most extreme ecological imbalance and highest functional activity of CRC drivers. But this subgroup exists equally in eoCRC and loCRC.

---

## 6. Keystone Species and Network Fragility

### Finding
F. prausnitzii is the top keystone species in both eoCRC (disruption score 0.124) and young controls (0.194). Removing it increases network fragmentation by 2 components in both groups.

P. stomatis is the #2 keystone in eoCRC (disruption 0.082) but is not in the top 5 for controls. An oral pathobiont becoming architecturally important to the disease network is biologically notable.

Removing butyrate producers from the healthy control network caused co-clusterings to appear (0→2 in yCtrl), suggesting that butyrate producers provide a structural buffer against oral-CRC integration. The eoCRC network was unchanged by butyrate producer removal — it already resembles the "buffering removed" state.

### Assessment: SUGGESTIVE
The keystone analysis is descriptive and the fragility simulation is based on one network topology. Interesting for the discussion but not a standalone finding.

---

## 7. Summary of What is Novel

1. **Dual-module framework for CRC dysbiosis** — framing CRC as a balance between two anticorrelated modules rather than a species list. Novel framing, robust across cohorts, not eoCRC-specific.

2. **Network diffusion analysis applied to CRC microbiome** — entirely new analytical approach for this field. Reveals directional influence flow invisible to abundance analysis. eoCRC shows extreme oral-to-CRC influence (ratio 40:1, permutation p<0.001), significantly higher than sample-size-matched loCRC (median 16.0 vs 5.4, Wilcoxon p=0.001, Cohen's d=0.46).

3. **eoCRC has more stereotyped microbial topology than loCRC** — eoCRC network topology is highly reproducible (CV=0.33 across leave-one-study-out), while loCRC at matched sample sizes shows high topological variability (CV=1.50). This suggests eoCRC converges on a uniform ecological architecture while loCRC encompasses diverse microbial ecologies.

4. **Gemella-P. stomatis co-occurrence replicates in 5/5 independent studies** — most microbiome network papers don't demonstrate within-study replication.

5. **P. stomatis as a keystone species in eoCRC** — an oral pathobiont becoming structurally important to the disease network.

6. **Compositional vs topological dysbiosis distinction** — eoCRC and loCRC are compositionally identical (same species, same abundances, same pairwise co-occurrence rates, same per-sample module scores) but topologically distinguishable (different network diffusion properties, different community structure). This concept — that identical species compositions can harbor different ecological architectures — is new to the CRC microbiome field.

---

## 8. Summary of What is NOT eoCRC-Specific

1. Species-level enrichment/depletion (same oral pathobionts, same butyrate depletion)
2. Dual-module ratio (same separation of CRC from controls in both age groups)
3. Per-sample co-occurrence scores (all p>0.4 for eoCRC vs loCRC)
4. Integrated patient proportion (25.6% vs 21.9%, p=0.47)
5. Functional pathway elevation in integrated patients (same pattern in loCRC)

---

## 9. Honest Limitations

1. **Sample size asymmetry:** eoCRC n=78 vs loCRC n=562. Network methods behave differently at different sample sizes. The full-cohort eoCRC-vs-loCRC diffusion contrast (40:1 vs 2.3:1) overstates the true difference. At matched sample sizes, the distinction remains significant (p=0.001) but is moderate (Cohen's d=0.46), with 17% of loCRC subsamples exceeding the eoCRC median.

2. **Prevalence asymmetry:** Oral pathobionts are at 1-10% prevalence in controls, limiting SpiecEasi's ability to estimate edges. Absence of co-clustering in controls may partly reflect inability to detect edges rather than true ecological absence.

3. **Cross-sectional design:** Cannot determine temporal direction. Tumor microenvironment may recruit oral bacteria rather than oral bacteria causing tumors.

4. **Public data limitations:** No tumor staging, MSI status, or tumor sidedness available for most samples. Cannot link network findings to clinical outcomes.

5. **P. stomatis-B. fragilis correlation does not replicate within YachidaS_2019** (rho=0.015), the largest individual cohort.

6. **Group-level vs per-sample:** The community detection and network diffusion findings are group-level network properties. Per-sample module scores (co-presence, co-abundance, integration ratio) do not distinguish eoCRC from loCRC (all p>0.4). The topological distinction exists at the network level but does not translate into measurable per-sample differences in species co-occurrence.

---

## 10. Files and Data Locations

### Saved R objects (in ~/eoCRC_analysis/results/networks/)
- `spieceasi_eoCRC.rds` — SpiecEasi model object, eoCRC
- `spieceasi_yCtrl.rds` — SpiecEasi model object, young controls
- `spieceasi_loCRC.rds` — SpiecEasi model object, loCRC
- `spieceasi_oCtrl.rds` — SpiecEasi model object, older controls
- `igraph_objects.rds` — igraph network objects for all 4 groups
- `species_annotations.rds` — oral_taxa, oral_origin, short_names lists
- `network_final_results.rds` — community detection results with stability
- `network_summary.csv` — summary statistics table
- `diffusion_results.rds` — all diffusion analysis results
- `matched_comparison.rds` — 100-iteration matched loCRC vs eoCRC LOSO distributions
- `module_scores.csv` — per-sample oral-CRC, butyrate, and dual-ratio scores

### Source data
- `~/eoCRC_analysis/species_filtered.rds` — 1282 x 220 species abundance matrix
- `~/eoCRC_analysis/all_samples_metadata.rds` — 1282 x 143 sample metadata
- `~/eoCRC_analysis/tse_pathways_full.rds` — 37945 x 1282 stratified pathway matrix

**NOTE:** Results may need to be re-saved if R session was interrupted. See save script (eoCRC_save_and_reproduce.R).
