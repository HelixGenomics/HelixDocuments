# Bayesian Evidence Scoring Pipeline

## Overview

Each genetic variant in a report is scored across 9 independent evidence databases. Evidence is combined using log-likelihood ratios, weighted by source quality and independence, producing a single evidence score per variant-disease pair.

This replaces the earlier approach of relying solely on PRS percentiles + ClinVar pathogenicity labels. The Bayesian pipeline allows weaker evidence from multiple sources to accumulate into meaningful signal, while strong single-source evidence (e.g., ClinVar 4-star pathogenic) can stand on its own.

## Evidence Sources

### 1. ClinVar (Star-Weighted)

ClinVar review status stars (0-4) weight the contribution of pathogenicity classifications:

| Stars | Meaning | Weight |
|:------|:--------|:-------|
| 4 | Practice guideline | 1.0 |
| 3 | Expert panel reviewed | 0.9 |
| 2 | Multiple submitters, no conflict | 0.7 |
| 1 | Single submitter | 0.4 |
| 0 | No assertion criteria | 0.1 |

A 4-star "pathogenic" classification contributes far more evidence than a 0-star one. This is enriched from the ClinVar VCF using our `build-clinvar-stars.py` tool.

### 2. CADD (Deleteriousness)

Combined Annotation Dependent Depletion score predicts variant deleteriousness. CADD scores > 20 (top 1% of variants) contribute positive evidence; scores < 10 contribute negative evidence. Phred-scaled scores are converted to log-likelihood ratios.

### 3. gnomAD (Population Frequency)

Rare variants (MAF < 0.01) in gnomAD contribute positive evidence for pathogenicity. Common variants (MAF > 0.05) contribute negative evidence. Uses ancestry-matched population frequencies when available.

### 4. AlphaMissense

DeepMind's protein structure model predicts missense variant impact. Scores > 0.564 (likely pathogenic threshold) contribute positive evidence.

### 5. ClinGen

Gene-disease validity classifications from the Clinical Genome Resource. "Definitive" and "Strong" evidence levels boost variant-disease associations in those genes.

### 6. DisGeNET

Disease-gene association scores from curated and text-mined sources. Gene Disease Association (GDA) scores > 0.3 contribute positive evidence.

### 7. GWAS Catalog

Genome-wide significant associations (p < 5e-8) from published GWAS studies. Effect sizes and sample sizes modulate evidence contribution.

### 8. PharmGKB

Pharmacogenomic evidence levels (1A through 4) for drug-gene interactions. Level 1A/1B associations are treated as strong evidence.

### 9. SNPedia

Curated variant annotations with magnitude scores (0-10). High-magnitude entries contribute supplementary evidence.

## Combination Method

Evidence from all sources is combined using the log-likelihood ratio framework:

```
LR_combined = LR_clinvar * LR_cadd * LR_gnomad * LR_alphamissense * ...
```

Each source produces an independent likelihood ratio (how much more likely is this variant pathogenic vs benign given this evidence). Ratios are multiplied across sources, then converted to a posterior probability using the prior prevalence of the condition.

### Independence Assumptions

Sources are treated as conditionally independent given the true pathogenicity status. This is approximately correct — ClinVar submissions may reference CADD or gnomAD, introducing some correlation, but in practice the independence assumption works well enough for evidence ranking.

### Carrier Detection

The pipeline identifies heterozygous carriers of autosomal recessive conditions. A carrier with one pathogenic allele in a recessive gene doesn't have the condition but is flagged for reproductive counseling. This is distinct from dominant conditions where a single pathogenic allele is sufficient.

## Implementation

The evidence pipeline is implemented in the server-side scoring engine (`server.js`). Evidence tables are pre-loaded from the unified SQLite database and queried per-variant during scoring. The combined evidence score is passed to specialist AI agents alongside PRS percentiles, enabling them to weigh strong single-variant evidence against polygenic background risk.

## Results

The Bayesian pipeline improved the quality of AI agent outputs in several ways:
- **Reduced false confidence**: Agents no longer report "elevated risk" for conditions with only PRS support and no variant evidence
- **Better carrier detection**: Heterozygous carriers of recessive conditions are correctly identified rather than being flagged as at-risk
- **Star-weighted ClinVar**: Agents can distinguish between well-reviewed pathogenic variants and single-submitter uncertain classifications
