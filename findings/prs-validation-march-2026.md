# PRS Validation Results — March 2026 (Batch v3)

## Summary

Validated PRS predictions against self-reported diagnoses from 29 PGP participants (batch v3 — 74 genomes scored, 29 with phenotype overlap). Used survey data with explicit Yes/No per condition to calculate proper AUC (Mann-Whitney U), sensitivity, and specificity.

**Key findings:**
- 12 conditions had sufficient data (n >= 3 per group) for validation
- Mean AUC: 0.540 across all validated conditions
- Cholesterol prediction stood out with AUC 0.770 (strong discrimination)
- 4 models flagged as performing worse than random (AUC < 0.50)

## Per-Condition Results

| Condition | PGS Model | AUC | n (has) | n (no) | Mean Has | Mean No | Sensitivity | Specificity | Flag |
|:----------|:----------|:----|:--------|:-------|:---------|:--------|:------------|:------------|:-----|
| Cholesterol | PGS001324 | 0.770 | 11 | 15 | 69.2 | 50.0 | 0.46 | 0.87 | GOOD |
| Atopic eczema | PGS001849 | 0.606 | 6 | 22 | 55.0 | 50.0 | 0.50 | 0.59 | |
| Colon polyps | PGS001846 | 0.600 | 8 | 20 | 60.9 | 50.0 | 0.62 | 0.60 | |
| Skin cancer (BCC) | PGS001323 | 0.583 | 6 | 22 | 53.0 | 50.0 | 0.50 | 0.64 | |
| Hypothyroidism | PGS001324 | 0.576 | 6 | 22 | 53.0 | 50.0 | 0.50 | 0.55 | |
| Cholelithiasis | PGS001174 | 0.547 | 3 | 25 | 56.9 | 50.0 | 1.00 | 0.52 | |
| Hemorrhoids | PGS001846 | 0.519 | 11 | 17 | 50.4 | 50.0 | 0.55 | 0.53 | |
| Migraine | PGS001282 | 0.503 | 9 | 19 | 51.5 | 50.0 | 0.44 | 0.74 | NO_BETTER_THAN_RANDOM |
| Inguinal hernia | PGS001854 | 0.493 | 3 | 25 | 55.9 | 50.0 | 0.33 | 0.64 | |
| Glaucoma | PGS001323 | 0.467 | 3 | 25 | 48.1 | 50.0 | 0.67 | 0.40 | |
| Asthma | PGS001849 | 0.409 | 7 | 22 | 42.6 | 50.0 | 0.29 | 0.64 | WORSE_THAN_RANDOM |
| Hypertension | PGS001838 | 0.406 | 8 | 20 | 44.4 | 50.0 | 0.38 | 0.55 | WORSE_THAN_RANDOM |

## Methodology

### Data Sources
- **Phenotypes**: PGP survey responses from 148 participants across 3 batches
  - Survey format 1: "Have you ever been diagnosed with...?" followed by comma-separated condition list
  - Survey format 2 (s23): Explicit Yes/No per condition (asthma, COPD, diabetes, etc.)
  - Profile traits: Conditions with diagnosis dates
- **PRS Scores**: Ensembled scores from batch v3 pipeline (post-blocklist, post-calibration)
- **Condition mapping**: 80+ self-reported conditions mapped to PRS trait names

### AUC Calculation
AUC computed using Mann-Whitney U statistic — equivalent to the probability that a randomly chosen affected individual scores higher than a randomly chosen unaffected individual. For traits where high_is_bad=false, scores are inverted before comparison.

### Limitations
- **Small sample sizes**: Most conditions have only 3-11 affected participants, making AUC estimates noisy
- **Self-report bias**: PGP conditions are self-reported, not clinically verified
- **Ascertainment bias**: PGP participants are not representative of the general population (skew toward educated, health-aware, primarily EUR ancestry)
- **Condition mapping**: Imperfect mapping between self-reported condition names and PRS trait names may introduce noise

## Blocklist Actions

Based on validation + manual distribution analysis:

| Model | Trait | Action | Reason |
|:------|:------|:-------|:-------|
| PGS001255 | Otosclerosis | **Blocked** | Mean percentile 16.3 across 108 reports (single model, severe downward bias) |
| PGS005162 | Colorectal cancer | **Blocked** | Mean percentile 17.3 across 99 reports (downward bias despite 4-model ensemble) |
| PGS001838 | Hypertension | Watch list | AUC 0.406 but small sample (n=8), 72 alternative models available |
| PGS001849 | Asthma | Watch list | AUC 0.409 but small sample (n=7), 67 alternative models available |
| PGS001323 | Glaucoma | Watch list | AUC 0.467 but tiny sample (n=3), 15 alternative models available |
| PGS001854 | Inguinal hernia | Watch list | AUC 0.493 but tiny sample (n=3), 6 alternative models available |

Watch list models will be reassessed when sample sizes reach n >= 20 per group.

## Validation Script

See [tools/validate-prs-phenotype.py](../tools/validate-prs-phenotype.py) for the full validation implementation.
