<p align="center">
  <img src="https://img.shields.io/badge/Tools-15_Scripts-8b5cf6?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/PGS_Models-2%2C826-06b6d4?style=for-the-badge" alt="PGS Models">
  <img src="https://img.shields.io/badge/Validated-PGP_%2B_OpenSNP-10b981?style=for-the-badge" alt="Validated">
  <img src="https://img.shields.io/badge/License-MIT-f59e0b?style=for-the-badge" alt="License">
</p>

<h1 align="center">🛠️ Research Tools</h1>

<p align="center">
  <strong>Working tools from the Helix pipeline — shared for reproducibility.</strong><br>
  <em>These aren't polished libraries. They're the actual scripts we use.</em>
</p>

<p align="center">
  <a href="https://helixsequencing.com">Platform</a> ·
  <a href="../journal/2026-03-prs-scoring-journey.md">Engineering Journal</a> ·
  <a href="../README.md">Back to Main</a>
</p>

---

## 🧬 Scoring & Ensemble

<table>
<tr>
<td width="240"><strong><code>prs-scorer-ensemble.py</code></strong></td>
<td>
PRS scorer with <strong>equal-weight trimmed mean</strong> ensemble averaging.<br><br>
<ul>
<li>Scores 2,826+ PGS Catalog models against imputed genotypes</li>
<li>EUR frequency-based ambiguous SNP resolution (A/T, C/G strand disambiguation)</li>
<li>Groups models by EFO trait, removes outliers (>1.0 from median z), simple average</li>
<li>Skips genome-wide models (>100K variants) — imputation bias inflates scores at scale</li>
<li>Skips <20 variant models (not true polygenic scores)</li>
<li>54 models blocklisted for non-discrimination, distribution bias, or phenotype mismatch</li>
<li>Frequency-based calibration with correction models</li>
</ul>
</td>
</tr>
<tr>
<td><strong><code>score-1kg-mega.py</code></strong></td>
<td>Scores 1000 Genomes Phase 3 samples against all PGS models to build reference population distributions.</td>
</tr>
<tr>
<td><strong><code>select-best-models.py</code></strong></td>
<td>Selects optimal PGS models per trait based on coverage, calibration, and validation performance.</td>
</tr>
<tr>
<td><strong><code>check-allele-alignment.py</code></strong></td>
<td>Checks allele alignment between scoring files and genotype data. Detects strand flips and ambiguous SNPs.</td>
</tr>
</table>

---

## ✅ Validation

<table>
<tr>
<td width="240"><strong><code>validate-pgp-chips.py</code></strong></td>
<td>
Validates PGP genotype files for chipset completeness and phenotype richness.<br><br>
<ul>
<li>Counts SNPs per file (full 23andMe = 574K–967K)</li>
<li>Identifies partial/broken files</li>
<li>Scores medical data richness (diagnoses, lab values, surveys)</li>
<li>Ranks participants by validation priority</li>
</ul>
</td>
</tr>
<tr>
<td><strong><code>validate-pgp.py</code></strong></td>
<td>Validates PRS predictions against actual PGP Harvard participant phenotypes. Compares genetic risk scores to self-reported diagnoses, lab values, and health outcomes.</td>
</tr>
<tr>
<td><strong><code>validate-opensnp-gpu-v2.py</code></strong></td>
<td>GPU-accelerated validation against OpenSNP phenotype data (~7,000 genomes).</td>
</tr>
</table>

---

## 📥 Data Collection

<table>
<tr>
<td width="240"><strong><code>deep-scrape-pgp.py</code></strong></td>
<td>Scrapes PGP Harvard participant profiles — surveys, phenotypes, health records, and genetic file metadata. Downloads top N participants ranked by data richness.</td>
</tr>
<tr>
<td><strong><code>fetch-pgs-metadata.py</code></strong></td>
<td>Downloads PGS Catalog metadata: trait names, EFO mappings, development methods, variant counts, and ancestry distributions.</td>
</tr>
<tr>
<td><strong><code>download-ancestry-pgs.py</code></strong></td>
<td>Downloads ancestry-specific PGS scoring files and population distributions from the PGS Catalog FTP.</td>
</tr>
<tr>
<td><strong><code>generate-trait-descriptions.py</code></strong></td>
<td>
Generates plain-English trait descriptions for PGS models using Claude.<br><br>
<ul>
<li>Produces layman-friendly descriptions (<40 words each)</li>
<li>Covers 1,721 unique traits</li>
<li>Batches 40 traits per API call, saves incrementally</li>
<li>Output used in consumer-facing PRS scorecards</li>
</ul>
</td>
</tr>
</table>

---

## 🔬 Analysis

<table>
<tr>
<td width="240"><strong><code>smart-condition-mapper.py</code></strong></td>
<td>Maps PGS trait names to clinical categories using NLP and keyword matching.</td>
</tr>
<tr>
<td><strong><code>build-pipeline-gpu.py</code></strong></td>
<td>GPU-accelerated PRS scoring pipeline for batch processing large cohorts.</td>
</tr>
<tr>
<td><strong><code>find-missing-ancestry.py</code></strong></td>
<td>Identifies PGS models missing ancestry-specific population distributions.</td>
</tr>
</table>

---

## 💡 Key Findings

Discoveries from building and validating the pipeline. These cost us real time and money to figure out — sharing so you don't have to.

### Genome-Wide Model Inflation

> Models with >100K variants produce **systematically inflated z-scores** on consumer chip data processed through Beagle imputation.

| Variant Range | Models | Avg \|z\| | Extreme (>2σ) |
|:---|:---:|:---:|:---:|
| 100 – 50K | 1,599 | 0.74 | ~5% ✅ |
| 100K – 500K | 105 | 1.56 | 28.6% ⚠️ |
| 500K – 1M | 437 | 1.55 | 28.8% ⚠️ |
| 1M+ | 614 | 2.50 | 45.8% ❌ |

**Root cause:** Tiny per-variant imputation bias compounds across hundreds of thousands of variants. A 1.27M-variant heart failure model scored 99.9th percentile for nearly everyone we tested.

**Fix:** Exclude >100K variant models from ensemble scoring entirely.

---

### ClinVar Indel False Positives

> Consumer DNA chips **cannot reliably call insertions/deletions**. Genotypes showing "I", "D", "II", "DD" produce false pathogenic calls.

We found 30+ MECP2 "Rett syndrome" pathogenic variants from a single 23andMe file — all with "I" or "D" genotypes. None were real.

**Fix:** Only accept `/^[ACGT]{1,2}$/` genotypes from consumer chip data.

---

### gnomAD Multi-Row Query Trap

> gnomAD has **multiple rows per rsID** (different transcript consequences). Using `LIMIT 1` can return the rare allele frequency instead of the common one.

Example: rs2815822 has rows with AF = 0.000007, 0.87, and 0.88. `LIMIT 1` returns 0.000007, making an 89%-frequency variant look rare — and leaving it classified as "Pathogenic" in ClinVar.

**Fix:** Use `MAX(af)` to get the highest frequency for any rsID.

---

### Equal-Weight Ensemble > Coverage-Weighted

> Coverage-weighted averaging lets genome-wide models (100% coverage) **dominate** over smaller, more targeted models that are actually more accurate.

Equal-weight trimmed mean — every model gets one vote, outliers removed — is the standard meta-analysis approach and proved more robust in our validation against PGP phenotypes.

---

<p align="center">
  <sub>Part of <a href="https://github.com/HelixGenomics/helix-open-research">Helix Open Research</a> · <a href="https://helixsequencing.com">helixsequencing.com</a></sub>
</p>
