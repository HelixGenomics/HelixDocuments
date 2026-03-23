<p align="center">
  <strong>🧬 HELIX SEQUENCING — ENGINEERING JOURNAL</strong><br>
  <em>Entry 1 · March 2026</em>
</p>

---

# How We Built Accurate Polygenic Risk Scores from a $99 DNA Kit

> **TL;DR:** We scored 2,800+ disease risk models from a consumer DNA file. The naive version was wildly wrong. Three weeks of "improvements" made it worse. Going back to basics — with honest confidence tiers — produced the best results.

---

## 🎯 The Goal

Take a 23andMe file. Score it against every published polygenic risk score in the [PGS Catalog](https://www.pgscatalog.org/) (3,550+ models). Tell someone what their DNA says about thousands of conditions. Make it accurate. Make it fast.

*It was harder than we expected.*

---

## 💥 What Went Wrong

Our first version produced garbage — users showed 99th percentile cancer risk everywhere. Scores were dramatic, bimodal, and **wrong**. Three compounding issues:

### 1. Coverage Mismatch
Consumer chips capture ~600K SNPs. Some PGS models need 1M+. We scored an 11%-coverage model the same as a 98%-coverage one.

> **Fix:** [Beagle 5.5](https://faculty.washington.edu/browning/beagle/) imputation expands ~600K → ~28M variants. Coverage jumps from ~30% to ~95%. Plus per-model confidence tiers so users know when to trust a score.

### 2. Ancestry Bias
One reference distribution for everyone. African ancestry scored against European values → systematically wrong percentiles.

> **Fix:** 14 ancestry-informative markers auto-detect population. Scores compared against matched distributions (EUR, AFR, EAS, SAS, AMR) from all 2,504 samples in 1000 Genomes Phase 3.

### 3. Strand-Ambiguous SNP Inflation
A/T and C/G pairs read in either orientation. Complement-matching inflated scores — Height PRS had 477/3290 ambiguous variants with 99.8% positive weights → **z-score of 14+**.

> **Fix:** Direct allele matching only. No complement flipping, ever. One change eliminated entire categories of false-positive extremes.

---

## 💸 The Expensive Mistake

Three weeks. ~$200 in GPU compute. We tried PRS-CSx (Bayesian LD-adjusted weights) and Ridge regression correction models on top.

**The correction models pulled everything toward the mean.**

| Before correction | After correction |
|:--:|:--:|
| 99th percentile (genuine high risk) | 55th percentile |
| 5th percentile (genuine protection) | 45th percentile |
| Clinically informative | Clinically useless |

Scores became flat. Featureless. A 99th percentile T1D risk driven by HLA variants (>50% of T1D genetic variance) got smoothed to "slightly above average." We had optimized for comfort over accuracy.

> **The lesson:** Extreme scores are not bugs. They're the whole point. The problem was showing low-confidence and high-confidence scores identically — not that scores were "too extreme."

We rolled everything back.

---

## 📊 Validation

Tested against **4,257 OpenSNP samples** and **84 Harvard Personal Genome Project genomes**:

| Trait | Metric | Our Result | Gold Standard (UK Biobank WGS) |
|-------|--------|:----------:|:------------------------------:|
| **Height** | Pearson r | 0.107 (p=0.044) | r ≈ 0.45–0.50 |
| **Red hair** | AUC | 0.67 | — |
| **Black hair** | AUC | 0.63 | — |
| **Eye colour** | AUC | 0.54 ⚠️ | ~0.95 (with rs12913832) |

Eye colour failure was instructive: missing **one SNP** (rs12913832, the strongest eye colour predictor) reduced classification to random chance. Coverage matters more than model sophistication.

---

## 💰 What It Cost to Build

| Item | Cost |
|------|-----:|
| GPU compute (Vast.ai) | $200–400 |
| Claude API (AI reports) | $100–300 |
| VPS hosting | $20/mo |
| Domain + SSL | $15/yr |
| **Total** | **~$400–800** |

The most expensive thing wasn't compute — it was three weeks of engineering time on models we threw away.

---

## ✅ Where We Landed

<table>
<tr><td width="50%">

**What works**
- 2,800+ PRS models with proper allele alignment
- Per-model confidence tiers
- Beagle 5.5 imputation (600K → 28M variants)
- 16 integrated databases (ClinVar, gnomAD, CPIC, AlphaMissense...)
- 6 parallel AI agents (sex-aware health narratives)
- <2 min standard / ~3-5 min with imputation
- Zero data retention + cryptographic deletion certs

</td><td>

**Still working on**
- Non-EUR ancestry accuracy
- Trio analysis for rare disease families
- Drug interaction deep-dives

</td></tr>
</table>

---

## 🔬 Update: Batch v3 & Bayesian Evidence Pipeline (Late March 2026)

Since the initial journal entry, significant progress:

### Bayesian Evidence Scoring
We moved from PRS-only scoring to a **9-source Bayesian evidence pipeline**. Each variant is now evaluated across ClinVar (with star-weighted classifications), CADD, gnomAD, AlphaMissense, ClinGen, DisGeNET, GWAS Catalog, PharmGKB, and SNPedia. Evidence accumulates via log-likelihood ratios — weak signals from multiple independent sources can build meaningful evidence.

The biggest win: **ClinVar star weighting**. A 4-star "practice guideline" pathogenic classification now contributes 10x more evidence than a 0-star single-submitter entry. This alone reduced false-positive pathogenic variant callouts by ~40%.

### Batch v3 Results
- **74 genomes scored** from PGP participants
- **1,434 training pairs** generated (42 unique participants x 8 agent types)
- **Proper AUC validation**: cholesterol prediction at AUC 0.770, mean across 12 conditions at 0.540
- **Blocklist grown to 71 models** (from 54), including newly identified biased models for otosclerosis and colorectal cancer
- **8 specialist agents** (added Body domain, split from general)

### Agent Domain Split
Reports are now generated by 8 specialist AI agents: Cardio, Cancer+Immune, Metabolic, Neuro, Body, Protocols, Health Index, and Summaries. Each agent receives the full evidence package (PRS + variants + pharmacogenomics) but focuses on its clinical domain.

### Health Index Recalibration
The composite health index needed recalibration after adding the evidence pipeline — pathogenic variant counts were skewing scores too aggressively. We built a [simulation tool](../tools/sim-health-index.py) to iterate on weights until the population distribution looked right (median ~70, genuine high-risk at 30-50).

### What's Next
- Process next ~75 genomes after quota reset (Mar 27)
- Reassess watch-listed PRS models (hypertension, asthma) when sample sizes reach n >= 20
- Begin fine-tuning experiments on the 1,434-pair training corpus

---

## 🔑 The Takeaway

> **More sophisticated models don't mean better predictions. Less is more.**

Every attempt to "fix" scores by adding complexity made them worse. The PGS Catalog models are peer-reviewed and validated in cohorts of millions. Our job is to make them accessible from a $99 DNA kit — not to second-guess them.

Score accurately. Present honestly. Let the data speak.

---

<p align="center">
  <a href="https://www.helixsequencing.com">helixsequencing.com</a> · Privacy-first DNA analysis · Zero data retention
</p>
