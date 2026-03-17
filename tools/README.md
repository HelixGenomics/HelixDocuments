<p align="center">
  <img src="https://img.shields.io/badge/Python-3.10+-3776ab?style=for-the-badge&logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/PyTorch-GPU_Accelerated-ee4c2c?style=for-the-badge&logo=pytorch&logoColor=white" alt="PyTorch">
  <img src="https://img.shields.io/badge/License-MIT-f59e0b?style=for-the-badge" alt="MIT">
</p>

<h1 align="center">🧬 Helix PRS Tools</h1>

<p align="center">
  <strong>Standalone Python scripts for polygenic risk score research.</strong><br>
  <em>All tools work with publicly available data. No proprietary datasets required.</em>
</p>

---

## ⭐ Headline Tool: PGP Deep Scraper

> **`deep-scrape-pgp.py`** — The only open-source tool that collects **real genotype data paired with real medical records** from a public dataset.

The [Harvard Personal Genome Project](https://my.pgp-hms.org/) has **943+ participants** who publicly share both their DNA data and their complete medical histories. This is the gold standard for PRS validation — real chip data, real diagnoses, real outcomes.

Our scraper extracts:
- **Structured disease diagnoses** mapped to ICD categories
- **Survey responses** (demographics, family history, medications)
- **Genotype file locations** for direct download
- **Phenotype JSON** ready for automated PRS validation

No other public dataset gives you genotypes + clinical outcomes at this scale. Most PRS validation uses self-reported traits (height, hair colour). This gives you **doctor-confirmed diagnoses** — cancer, cardiovascular disease, autoimmune conditions, metabolic disorders.

```bash
python3 deep-scrape-pgp.py --output pgp-phenotypes.json
# → 943 profiles, 579 with medical conditions, 911 unique diagnoses
```

---

## Prerequisites

- Python 3.10+
- Common: `numpy`, `scipy`, `requests`
- GPU scripts: `torch` with CUDA (tested on RTX 3060 Ti, 8GB VRAM)
- [PGS Catalog](https://www.pgscatalog.org/) API access (free)
- [OpenSNP Frei 2024](https://zenodo.org/records/10715132) dataset (4,257 QC'd genomes) for validation
- [Harvard PGP](https://my.pgp-hms.org/) profiles for phenotype validation

---

## 📦 Data Collection

| Script | What It Does |
|:-------|:-------------|
| **`deep-scrape-pgp.py`** | Scrape 943+ PGP profiles — genotypes + medical records + surveys. Extracts structured disease diagnoses for PRS validation. **The key tool for clinical ground truth.** |
| `download-ancestry-pgs.py` | Bulk download PGS Catalog scoring files. Converts to compact gzipped TSV. 20 parallel workers. |
| `fetch-pgs-metadata.py` | Fetch ancestry metadata for all PGS models via REST API. Builds per-trait, per-ancestry model selection map. |
| `find-missing-ancestry.py` | Find multi-ancestry models not yet downloaded. Identifies models with >5% non-EUR representation. |

## 🔬 Mapping & QC

| Script | What It Does |
|:-------|:-------------|
| **`smart-condition-mapper.py`** | Map free-text clinical conditions to PRS trait names. Handles abbreviations (CAD, HTN, COPD), qualifier stripping, keyword matching against 170+ patterns. Built-in test suite. |
| `check-allele-alignment.py` | Verify allele alignment between PGS effect alleles and genotype data. Catches strand flips and complement mismatches that cause systematic scoring errors. |

## ⚡ Scoring

| Script | What It Does |
|:-------|:-------------|
| `score-1kg-mega.py` | Score all 2,504 1000 Genomes samples across 5 ancestries (EUR/AFR/EAS/SAS/AMR). Builds population distributions for percentile estimation. |
| **`build-pipeline-gpu.py`** | Full GPU pipeline — VCF parsing, allele alignment, batched scoring against 3,550 models (2.37B weights). Built for cloud GPU instances. |
| `validate-opensnp-gpu-v2.py` | GPU validation against 4,257 OpenSNP genomes. 6.2GB sparse weight matrix, chunked sparse multiplication. Reports Pearson r, AUC, per-ancestry distributions. |

## ✅ Validation & Selection

| Script | What It Does |
|:-------|:-------------|
| **`validate-pgp.py`** | Validate PRS against PGP phenotypes — AUC for disease case/control, Pearson r for continuous traits. Uses the medical records from `deep-scrape-pgp.py`. |
| `select-best-models.py` | Combine PGP + OpenSNP results → select best PGS model per trait per ancestry → `best-models.json` for production. |

---

## 🔄 Data Flow

```
PGS Catalog ──[download-ancestry-pgs.py]──► Compact scoring files
     │
     └──────[fetch-pgs-metadata.py]───────► Per-ancestry model map

Harvard PGP ──[deep-scrape-pgp.py]────────► Phenotype JSON (943 profiles)
     │                                        │
     │                                        ├── 579 with medical conditions
     │                                        ├── 911 unique diagnoses
     │                                        └── Genotype download links
     │
     └──────[smart-condition-mapper.py]───► Standardized trait mappings

1000 Genomes ─[score-1kg-mega.py]─────────► Population distributions (5 ancestries)

OpenSNP ──────[validate-opensnp-gpu-v2.py]► Validation metrics (r, AUC)

PGP Genomes ──[validate-pgp.py]───────────► Disease AUC per model

All results ──[select-best-models.py]─────► Production best-models.json
```

---

## ⚠️ Notes

- Scripts use hardcoded paths from our build environment — update for your setup
- GPU scripts require PyTorch + CUDA (tested on RTX 3060 Ti, 8GB VRAM)
- The mega SQLite database (`pgs-mega.db`, 174GB) is not included — build from PGS Catalog bulk downloads
- All tools use **publicly available datasets only**

---

## 📄 License

MIT — use freely for any purpose.

<p align="center">
  <a href="https://github.com/HelixGenomics/helix-open-research">← Back to Helix Open Research</a>
</p>
