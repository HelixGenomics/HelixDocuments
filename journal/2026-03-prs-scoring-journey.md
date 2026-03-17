# Building Accurate Polygenic Risk Scores from Consumer DNA Data: An Engineering Journal

**Entry 1 — March 2026**
*Helix Sequencing Research & Development*

---

## The Problem We Set Out to Solve

When we started Helix Sequencing, the goal was straightforward: take a consumer genotype file from 23andMe or AncestryDNA (~600K–950K SNPs), score it against every relevant polygenic risk score model in the PGS Catalog, and tell the user something meaningful about their genetic risk for thousands of conditions.

The PGS Catalog contains 3,550+ published models covering everything from Type 2 Diabetes to height to hair colour. Each model is a set of genetic variants with effect weights — multiply your genotype by the weight, sum them up, and you get a raw polygenic score.

Simple in theory. In practice, the journey from raw score to accurate percentile has been one of the most technically demanding engineering problems we've encountered — and the most expensive.

---

## Phase 1: The Naive Approach (January 2026)

### What we built

Our first implementation was textbook: for each of the 3,550 PGS models, iterate through every variant, look up the user's genotype by rsID, multiply by the effect weight, and sum. The raw score gets converted to a percentile by comparing against a reference distribution.

```
raw_score = SUM(user_dosage × effect_weight)
percentile = binary_search(raw_score, sorted_reference_scores)
```

The reference distribution came from scoring all 2,504 samples in the 1000 Genomes Phase 3 dataset. We built distributions for 1,171 of the 3,550 models (those where at least 10 reference samples had enough variant coverage).

### What we got

Scores that looked dramatic. Many users showed 99th+ percentile for multiple cancer risks, or <1st percentile for common conditions. The distributions were heavily bimodal — scores clustered at the extremes.

### What was wrong

Several compounding issues:

1. **No ancestry matching.** We used a single EUR-based reference distribution for all users. A user with African ancestry scored against European reference scores will get systematically wrong percentiles because linkage disequilibrium patterns differ across populations.

2. **Allele alignment errors.** Some PGS models use the alternate allele as the effect allele, others use the reference allele. Without careful strand alignment, a proportion of variants were scored backwards — adding risk weight where it should subtract, and vice versa.

3. **Coverage gaps.** Consumer chips capture ~600K–950K SNPs. Large PGS models can have over 1 million variants. Coverage ranged from 1% to 99% depending on the model. We were scoring models with 11% coverage alongside models with 98% coverage, presenting both with equal confidence.

4. **No z-score fallback.** Models without reference distributions (2,379 of 3,550) used a crude z-score approximation that assumed normality — which polygenic scores rarely are.

---

## Phase 2: Imputation + Expanded Scoring (February 2026)

### The hypothesis

If low variant coverage is the problem, expand coverage through genotype imputation. Beagle 5.5 can impute ~600K chip SNPs up to ~28 million variants using the 1000 Genomes Phase 3 reference panel.

### What we built

A full imputation pipeline on a Vast.ai build server (256 vCPUs, 503GB RAM, RTX 3060 Ti):

1. Convert consumer chip data to VCF format
2. Run Beagle 5.5 per-chromosome imputation against 1000G reference
3. Score the imputed genotypes against all 3,550 models
4. Build per-ancestry reference distributions (EUR, AFR, EAS, SAS, AMR)

We also built a 174GB SQLite database (`pgs-mega.db`) containing all 2.375 billion variant weights across 3,550 models, indexed by both rsID and chromosomal position.

### The scale of it

Scoring 4,257 OpenSNP samples × 3,550 models required a 6.2GB sparse weight matrix (1.2M variants × 3,550 models, 1.24 billion non-zero entries). GPU-accelerated scoring on the RTX 3060 Ti processed the entire matrix in 7.4 minutes using chunked sparse CSR tensor multiplication via PyTorch.

### Validation results

We validated against OpenSNP phenotype data (4,257 QC'd samples from the Frei 2024 processed dataset):

**Height** (the gold standard PRS validation trait):
- Best model: PGS003888
- Pearson r = 0.107 (EUR), 0.095 (all ancestries)
- r² = 0.0115
- Coverage: 98.6% (1,053,219 matched variants)

For context, published height PRS studies in controlled cohorts like UK Biobank achieve r² ≈ 0.20–0.25 with whole-genome sequencing data. Our r² of 0.01 from consumer chip data with imputation was disappointing, but the correlation was statistically significant (p = 0.044).

**Hair colour** (binary classification):
- Best AUC: 0.67 for red hair prediction (PGS001094)
- Black hair prediction: AUC 0.63 (PGS001096)
- These are decent for consumer chip data — hair colour is highly polygenic with a few large-effect loci

**Eye colour** (using IrisPlex SNPs):
- Only 4 of 6 IrisPlex SNPs found in the genotype data
- AUC for blue vs brown: 0.54 — barely above random
- Missing rs12913832 (HERC2/OCA2), the single strongest eye colour predictor, crippled this entirely

---

## Phase 3: The PRS-CSx Experiment (March 2026)

### The hypothesis

Raw PGS Catalog weights are typically GWAS beta coefficients — they don't account for linkage disequilibrium, where nearby correlated SNPs double-count the same genetic signal. Bayesian methods like PRS-CSx can produce posterior effect sizes that are theoretically 10–30% more accurate by jointly modelling GWAS data from multiple populations.

### What we did

We ran PRS-CSx on 33 quantitative traits (blood biomarkers, BMI, height, blood pressure) using UK Biobank GWAS summary statistics and 1000 Genomes LD reference panels. Each model produced ~1.1 million posterior effect sizes across 5 ancestry groups.

This ran on our Vast.ai build server over multiple days.

### What we built around it

A complete GPU-accelerated pipeline with multiple iterations:

- **Pipeline v1:** Sequential scoring, 64 parallel workers. Hit memory limits — BrokenPipeError when trying to hold 2.37 billion variants in memory across worker processes.
- **Pipeline v2:** Parallel loading with 64 workers. Same memory crash at ~98 million variants loaded.
- **Pipeline v3:** Batched scoring — split 3,550 models into 17 batches of ~150M variants each. This worked. Scored 84 PGP VCFs × 3,550 models successfully.

### What went wrong

The PRS-CSx posterior weights didn't improve our validation metrics over raw PGS Catalog weights. In some cases, they made things worse.

The fundamental issue: PRS-CSx is designed for large biobank-scale GWAS summary statistics with millions of participants. Our validation cohort was ~4,257 OpenSNP samples and ~84 PGP genomes. The Bayesian shrinkage was tuned for a scale of data we simply didn't have for validation.

We spent approximately 2–3 weeks of continuous GPU compute on this experiment across multiple Vast.ai instances.

---

## Phase 4: The PGP Validation Push (March 2026)

### The hypothesis

More validation data with real medical records would let us identify which PGS models actually predict disease. The Harvard Personal Genome Project (PGP) has 943 participants with publicly shared genotypes AND medical diagnoses — the ideal validation dataset.

### What we built

We scraped all 943 PGP profiles, extracting:
- 579 participants with medical conditions (61%)
- 911 unique medical conditions
- 6,093 total condition entries

We mapped these to PGS models and planned to impute + score all of them.

### What actually happened

Of 569 PGP genomes we attempted to impute, **538 failed** (95% failure rate). The dominant failure mode was bcftools merge/concat errors and Beagle memory exhaustion. Only 31 genomes imputed successfully in a 9.8-hour batch run.

The failures were mostly due to:
- Heterogeneous genotype formats across PGP participants (different chip versions, different coordinate systems)
- Memory constraints on imputation (Beagle with -Xmx2g was insufficient for some genomes)
- bcftools indexing timeouts on the temporary VCF files

We ultimately got ~84 PGP genomes scored across all 3,550 models. With only ~84 samples, our validation power for any individual disease was extremely limited — most conditions had fewer than 20 cases.

---

## Phase 5: The Correction Model Trap (March 2026)

### What we tried

Seeing that raw percentiles were clustering at extremes (many 99+ and <1 scores), we attempted to train Ridge regression correction models to "recalibrate" the scores. The idea was to use features like:

- Raw PRS score
- Variant coverage percentage
- Ancestry proportions
- Model size (number of variants)
- Weight balance (positive vs negative)

And train against empirical percentiles from the PGP/OpenSNP validation data.

### The critical mistake

**The correction model pulled everything toward the mean.**

When you train a regression to predict percentiles from raw scores using a small, noisy validation dataset, the optimal strategy for the model is to be conservative — predict everyone as roughly average. This minimizes error on the training data but destroys the clinical utility of the scores.

The result: conditions where a user genuinely had elevated genetic risk showed as 50th–60th percentile. Conditions where they had genuinely low risk also showed as 40th–50th percentile. The scores became uninformative — a flat, featureless landscape where nothing stood out.

### The realization

**Extreme scores are not necessarily wrong.** A user showing 99th percentile for Type 1 Diabetes risk on a well-validated model (PGS with AUC > 0.80 for T1D, driven by HLA-DR/DQ variants that account for >50% of genetic variance) is probably actually at the 99th percentile.

The problem was never that the percentiles were too extreme. The problem was that we were presenting low-confidence scores (from models with 11% coverage or no validation data) with the same weight as high-confidence scores (98% coverage, strong AUC).

---

## What We Spent

### GPU Compute (Vast.ai)

Our primary build server was a 256 vCPU / 503GB RAM / RTX 3060 Ti instance. At typical Vast.ai rates:

| Period | Usage | Approximate Cost |
|--------|-------|-----------------|
| Initial pipeline development (Feb 2026) | ~1 week continuous | ~$50–80 |
| PRS-CSx weight computation (Feb–Mar 2026) | ~2 weeks continuous | ~$100–200 |
| PGP batch imputation attempts (Mar 2026) | ~10 hours | ~$5–10 |
| OpenSNP GPU validation runs (Mar 2026) | ~33 minutes per run × multiple | ~$5–10 |
| RTX 5070 secondary server (Mar 2026) | ~$0.087/hr, intermittent | ~$20–40 |
| **Estimated total GPU compute** | | **$200–400** |

### API and Data Costs

| Item | Cost |
|------|------|
| Anthropic Claude API (AI report generation) | ~$100–300 |
| VPS hosting (production server) | ~$20/month |
| Domain + SSL | ~$15/year |
| Coinbase Commerce integration | Transaction fees only |
| **Total operational** | **~$200–500** |

### Crypto Payments Received (for context)

Recent transaction history shows incremental revenue of approximately $355 in crypto payments between late February and mid-March 2026, across Coinbase Polygon and Base networks.

### The expensive lesson

The most costly experiment wasn't the GPU time — it was the 2–3 weeks spent on the PRS-CSx correction pipeline that we ultimately rolled back. The compute cost was modest (~$200), but the opportunity cost of engineering time was significant.

---

## Where We Are Now (March 2026)

After all of this, we made the decision to **go back to basics**. Here's our current position:

### What works

1. **3,550 PGS Catalog models scored correctly** with proper allele alignment
2. **Ancestry detection** with per-population model selection (EUR, AFR, EAS, SAS, AMR)
3. **Beagle 5.5 imputation** from ~700K to ~28M variants (for individual users, not batch)
4. **Population reference distributions** from 1000 Genomes (2,504 samples, 5 populations)
5. **Raw percentiles that reflect genuine genetic variation** — extreme scores are preserved when the underlying model is well-powered
6. **ClinVar scanning** (400K+ pathogenic variants), **pharmacogenomics** (34 CPIC genes), **SNPedia annotations** (16K+ variants)
7. **AI-powered longevity protocols** from 6 parallel specialized Claude instances
8. **Zero data retention architecture** — all user data deleted after report generation

### What doesn't work yet

1. **Low-coverage model confidence scoring.** We present models with 11% coverage the same as models with 98% coverage. Users need confidence indicators.
2. **Non-EUR ancestry accuracy.** Most PGS models were developed in European cohorts. Prediction accuracy drops significantly for non-EUR users, even with ancestry-matched distributions.
3. **Phenotype validation at scale.** With only ~84 scored PGP genomes, we can't definitively say which models perform best for which traits.
4. **The eye colour problem.** Missing key SNPs (rs12913832) from some chip versions makes certain predictions impossible regardless of model quality.

### The philosophical shift

We stopped trying to make every score look moderate and clinically comfortable. The PGS Catalog models are published, peer-reviewed, and validated in large biobanks. Our job is to score them accurately and present them honestly — not to smooth them into a bell curve that feels reassuring but isn't informative.

A 99th percentile Type 2 Diabetes PRS from a model validated at AUC 0.68 in UK Biobank is genuinely meaningful information. A 15th percentile height PRS that corresponds to being shorter than average is validated by the r = 0.107 correlation we measured. These are real signals, not artifacts.

---

## The Roadmap

### Near-term (Q2 2026)

1. **Per-model confidence tiers.** Tag each scored trait as high/moderate/low/unreliable confidence based on variant coverage, validation AUC, and model size. Display this prominently in reports.
2. **Expand PGP validation.** Fix the batch imputation pipeline (increase Beagle memory allocation, handle format heterogeneity better) and score all 943 PGP genomes.
3. **OpenSNP + PGP combined validation.** Use the full ~5,000 sample validation cohort to select the best model per trait per ancestry.

### Medium-term (Q3–Q4 2026)

4. **Bayesian weight recomputation** using PRS-CSx/LDpred2 — but only for traits where we have enough validation data to confirm improvement. No more speculative bulk recomputation.
5. **RICE integration** for combined rare + common variant scoring (25.7% improvement reported in literature).
6. **GenoPred pipeline adoption** for standardized multi-method PRS comparison.

### Long-term

7. **GPU cluster scaling** with 20+ Vast.ai instances for distributed PRS-CSx computation across all validated traits.
8. **UKBB/All of Us access** for gold-standard validation (requires institutional application).
9. **Trio analysis capability** for rare disease family studies — directly relevant to the Trisomy 9 research that motivated the entire project.

---

## Lessons Learned

1. **More data does not automatically mean better scores.** Adding PRS-CSx posterior weights, Ridge regression corrections, and ensemble methods can actively degrade accuracy when your validation cohort is small.

2. **Extreme percentiles are features, not bugs.** Polygenic scores are supposed to spread people across a distribution. Clamping at 99+ or <1 is fine if the underlying model is well-powered. The mistake is presenting all models with equal confidence.

3. **The validation bottleneck is the real constraint.** We have 3,550 models and 2.37 billion variant weights. What we lack is ground truth — large numbers of genotyped individuals with confirmed phenotypes. Everything else (GPU compute, Bayesian methods, imputation) is engineering. The validation data is science.

4. **Imputation is fragile at scale.** A 95% failure rate on batch PGP imputation taught us that heterogeneous consumer chip data needs extensive per-sample QC before imputation, not just a batch pipeline.

5. **Start with what you can validate, then expand.** Instead of scoring all 3,550 models and hoping they're accurate, focus on the ~55 traits where PGP provides ground truth, validate those rigorously, and use that as the foundation.

---

*This is the first entry in an ongoing engineering journal. We'll publish updates as our validation pipeline matures and accuracy metrics improve.*

*All source code, pipeline scripts, and validation data referenced in this entry are available in our [GitHub repository](https://github.com/HelixGenomics/HelixSequencing).*

---

**Helix Sequencing** — Privacy-first DNA analysis. Zero data retention.
