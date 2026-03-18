# Research Tools

Scripts used in Helix Sequencing research. These are working tools from our pipeline — shared for reproducibility and collaboration.

## Scoring & Validation

### prs-scorer-ensemble.py
PRS scorer with equal-weight trimmed mean ensemble averaging. Key features:
- Scores 2,826 PGS Catalog models against imputed genotypes
- EUR frequency-based ambiguous SNP resolution (A/T, C/G strand disambiguation)
- Equal-weight trimmed mean ensemble: groups models by EFO trait, removes outliers (>1.0 from median z), simple average
- Skips genome-wide models (>100K variants) — we found systematic imputation bias inflates scores at this scale
- Skips <20 variant models from ensemble (not true polygenic scores)
- Frequency-based calibration with correction models

### validate-pgp-chips.py
Validates PGP genotype files for chipset completeness and phenotype richness. Checks:
- SNP count per file (full 23andMe chip = 574K-967K)
- Identifies partial/broken files
- Scores medical data richness (diagnoses, lab values, surveys)
- Ranks participants for validation priority

### validate-pgp.py
Validates PRS predictions against actual PGP participant phenotypes.

### validate-opensnp-gpu-v2.py
GPU-accelerated validation against OpenSNP phenotype data.

## Data Collection

### deep-scrape-pgp.py
Scrapes PGP Harvard participant profiles, surveys, phenotypes, and genetic file metadata. Downloads top N participants ranked by data richness.

### fetch-pgs-metadata.py
Downloads PGS Catalog metadata including trait names, EFO mappings, development methods, and variant counts.

### download-ancestry-pgs.py
Downloads ancestry-specific PGS scoring files and population distributions.

### generate-trait-descriptions.py
Generates plain-English trait descriptions for PGS models using Claude. Produces layman-friendly (<40 word) explanations for 1,700+ traits. Output used in consumer-facing PRS scorecards.

## Analysis

### score-1kg-mega.py
Scores 1000 Genomes Phase 3 samples against all PGS models to build reference population distributions.

### select-best-models.py
Selects optimal PGS models per trait based on coverage, calibration, and validation performance.

### check-allele-alignment.py
Checks allele alignment between scoring files and genotype data, detects strand flips.

### smart-condition-mapper.py
Maps PGS trait names to clinical categories using NLP and keyword matching.

### build-pipeline-gpu.py
GPU-accelerated PRS scoring pipeline for batch processing.

### find-missing-ancestry.py
Identifies PGS models missing ancestry-specific distributions.

## Key Findings

### Genome-Wide Model Inflation
Models with >100K variants produce systematically inflated z-scores on consumer chip data processed through Beagle imputation. Data analysis across 2,826 models showed:
- 100-50K variants: avg |z|=0.74, ~5% extreme (expected)
- 100K-500K variants: avg |z|=1.56, 28.6% extreme
- 1M+ variants: avg |z|=2.50, 45.8% extreme

Root cause: tiny per-variant imputation bias compounds across hundreds of thousands of variants. Fix: exclude >100K variant models from ensemble scoring.

### ClinVar False Positives on Consumer Chips
Consumer DNA chips cannot reliably call insertions/deletions. Genotypes showing "I", "D", "II", "DD" should be filtered — they produce false pathogenic calls (e.g., 30+ MECP2 Rett syndrome variants from a single chip file). Fix: only accept /^[ACGT]{1,2}$/ genotypes.

### gnomAD Multi-Row rsID Issue
gnomAD has multiple rows per rsID (different transcript consequences). Using `LIMIT 1` can return the rare allele frequency instead of the common one. Fix: use `MAX(af)` to get the highest frequency for any rsID, ensuring common variants are properly identified.

### Equal-Weight Ensemble Beats Coverage-Weighted
Coverage-weighted averaging lets genome-wide models (100% coverage) dominate over smaller, more accurate targeted models. Equal-weight trimmed mean gives every model one vote, with outlier removal. Standard meta-analysis approach, more robust to individual model bias.
