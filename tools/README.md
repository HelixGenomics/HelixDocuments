# Helix PRS Tools

Standalone Python scripts for polygenic risk score (PRS) research using publicly available data. These tools were developed during the Helix Sequencing PRS accuracy improvement pipeline and are shared for the benefit of the genomics research community.

## Prerequisites

- Python 3.10+
- Dependencies vary per script (see headers) — common: `numpy`, `scipy`, `torch` (for GPU scripts)
- Access to [PGS Catalog](https://www.pgscatalog.org/) API
- For validation: [OpenSNP Frei 2024 dataset](https://zenodo.org/records/10715132) (4,257 QC'd genomes)
- For phenotype validation: [Harvard PGP](https://my.pgp-hms.org/) profiles

## Tools

### Data Collection

| Script | Description |
|--------|-------------|
| `deep-scrape-pgp.py` | Scrape PGP profiles for medical records, survey responses, and phenotype data. Extracts structured disease diagnoses from 943+ participants for PRS validation. |
| `download-ancestry-pgs.py` | Bulk download PGS scoring files from the PGS Catalog FTP. Converts to compact gzipped TSV format (rsid/ea/oa/weight). 20 parallel workers. |
| `fetch-pgs-metadata.py` | Fetch ancestry metadata for all PGS models via the PGS Catalog REST API. Builds a per-trait, per-ancestry model selection map scoring each model on variant count, GWAS ancestry representation, and sample size. |
| `find-missing-ancestry.py` | Scan the PGS Catalog for multi-ancestry models not yet downloaded. Identifies models with >5% non-EUR ancestry representation. |

### Mapping & QC

| Script | Description |
|--------|-------------|
| `smart-condition-mapper.py` | Map free-text clinical conditions (from PGP profiles) to standardized PRS trait names. Handles abbreviation expansion (CAD, HTN, COPD), qualifier stripping, keyword matching against 170+ patterns, and multi-level confidence scoring. Includes a built-in test suite. |
| `check-allele-alignment.py` | Verify allele alignment between PGS effect alleles and genotype data. Catches strand flips and complement mismatches that cause systematic scoring errors. |

### Scoring

| Script | Description |
|--------|-------------|
| `score-1kg-mega.py` | Score all 2,504 1000 Genomes Phase 3 samples against PGS models from a mega SQLite database. Produces per-ancestry (EUR/AFR/EAS/SAS/AMR) population distributions for percentile estimation. Parallelized across chromosomes and models. |
| `build-pipeline-gpu.py` | Full GPU-accelerated PRS scoring pipeline. Handles VCF parsing, allele alignment, batched scoring against 3,550 models (2.37B variant weights), and distribution building. Designed for Vast.ai instances with 256+ vCPUs and GPU. |
| `validate-opensnp-gpu-v2.py` | GPU-accelerated PRS validation against 4,257 OpenSNP genomes. Loads a 6.2GB sparse weight matrix, performs chunked sparse matrix multiplication on GPU, validates against height/hair/eye phenotypes. Reports Pearson r, AUC, and per-ancestry distributions. |

### Validation & Selection

| Script | Description |
|--------|-------------|
| `validate-pgp.py` | Validate PRS scores against PGP phenotype data. Computes AUC for binary traits (disease case/control) and Pearson r for continuous traits (height, BMI). Maps PGP medical conditions to PGS models. |
| `select-best-models.py` | Combine PGP + OpenSNP validation results and select the best PGS model per trait per ancestry. Outputs `best-models.json` for production deployment. |

## Data Flow

```
PGS Catalog ──[download-ancestry-pgs.py]──> Compact scoring files
     │
     └──[fetch-pgs-metadata.py]──> Per-ancestry model selection map

Harvard PGP ──[deep-scrape-pgp.py]──> Phenotype JSON
     │
     └──[smart-condition-mapper.py]──> Standardized trait mappings

1000 Genomes ──[score-1kg-mega.py]──> Population distributions (5 ancestries)

OpenSNP ──[validate-opensnp-gpu-v2.py]──> Validation metrics (r, AUC)

PGP Genomes ──[validate-pgp.py]──> Disease AUC per model

All results ──[select-best-models.py]──> Production best-models.json
```

## Notes

- All scripts use hardcoded paths from our build environment. You'll need to update paths for your setup.
- GPU scripts require PyTorch with CUDA support and were tested on RTX 3060 Ti (8GB VRAM).
- The mega SQLite database (`pgs-mega.db`, 174GB) is not included — build it using the PGS Catalog bulk downloads.
- These tools work with publicly available datasets only (PGS Catalog, 1000 Genomes, OpenSNP, Harvard PGP).

## License

MIT — use freely for research purposes.
