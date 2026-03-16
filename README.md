# Helix Sequencing

Privacy-first DNA analysis platform. Upload a consumer genetic data file, get 960+ polygenic risk scores, pathogenic variant scanning, pharmacogenomics with star allele calling, and AI-powered longevity protocols — with zero data retention.

## Key Numbers

| Metric | Value |
|--------|-------|
| Polygenic Risk Scores | 3,500+ (PGS Catalog) |
| Curated PRS Traits | 44 (coverage-optimized, ML-corrected) |
| Ancestry Populations | 5 (EUR, AFR, EAS, SAS, AMR) + blended |
| ClinVar variants scanned | 400,000+ |
| SNPedia annotations | 16,000+ |
| Pharmacogenes (CPIC) | 34 genes, 110K+ diplotypes |
| Imputed variants | ~28M (Beagle 5.5) |
| Reference allele freqs | 3.5M variants × 6 populations |
| Analysis time | <2 min standard, 15-30 min imputed |
| Data retention | Zero |

---

## Polygenic Risk Scoring (PRS) — Deep Dive

### Overview

The PRS engine scores an individual's genetic predisposition across 960+ traits using published polygenic score models from the [PGS Catalog](https://www.pgscatalog.org/). A curated subset of 44 traits is presented in reports, selected for clinical relevance and chip/imputation coverage.

### Scoring Pipeline

#### 1. Variant Loading

Two genotype sources are supported:

- **Chip-only (~608K variants):** Direct genotype calls from consumer DNA chips (23andMe, MyHeritage, AncestryDNA, etc.). Parsed as hard-call alleles (e.g., `AG`, `TT`).
- **Imputed (~28M variants):** Beagle 5.5 imputation against 1000 Genomes Phase 3 reference panel. Stored in SQLite with continuous dosage values (0.0–2.0) and genotype posterior probability (`gp_max`).

```
Chip:     rsid → "AG"           (discrete allele pair)
Imputed:  rsid → {ref: "A", alt: "G", dosage: 1.03, gp_max: 0.94}
```

#### 2. Dosage Calculation

For each PGS variant, the effect allele dosage is computed:

- **Chip mode:** Count occurrences of the effect allele in the genotype string (0, 1, or 2).
- **Imputed mode:** Map the dosage field based on whether the effect allele matches `alt` (use dosage directly) or `ref` (use `2 - dosage`).

**Strand-ambiguous SNP handling (A/T and C/G pairs):**
- Only **direct allele matching** is used — no complement flipping.
- This prevents the systematic bias where complement matching would assign dosage=2 to all strand-ambiguous variants, inflating scores for PGS models with many such variants.
- All three layers (getDosage, getEurFreq, frequency file builder) enforce this consistently.

#### 3. Raw Score Computation

```
raw_score = SUM(dosage_i * weight_i)  for all matched variants
```

Where `weight_i` comes from the PGS Catalog effect weights.

#### 4. Percentile Calculation (Three-Tier)

**Primary — Empirical distribution (preferred):**
- Raw score is ranked against pre-computed scores from 503 EUR samples (1000 Genomes Phase 3).
- Population distributions are available for 6 groups: EUR, AFR, EAS, SAS, AMR, ALL.
- Binary search gives the exact percentile position.
- Clamped to [0.1%, 99.9%].

**Fallback — Frequency-based z-score:**
- When no population distribution exists, the expected mean and variance are computed from European allele frequencies:
  ```
  E[score] = SUM(2 * p_i * w_i)
  Var[score] = SUM(2 * p_i * (1-p_i) * w_i^2)
  ```
- Scale factor `matched/freqMatched` compensates for variants missing frequency data.
- Z-score clamped to [-3.5, +3.5], converted to percentile via normal CDF.

**Last resort — Population statistics:**
- Uses pre-computed population mean/std when available.

#### 5. Confidence Assessment

Each score receives a confidence rating based on variant coverage:

| Confidence | Chip Coverage | Imputed Coverage |
|-----------|--------------|-----------------|
| High | >= 70% | >= 60% effective |
| Moderate | >= 40% | >= 40% effective |
| Low | >= 15% | >= 20% effective |
| Unreliable | < 15% | < 20% effective |

**Imputation quality weighting:**
- High quality (GP >= 0.9): counts fully
- Medium quality (GP 0.7–0.9): counts at 50%
- Low quality (GP < 0.7): discounted
- Effective coverage = `(HQ + MQ*0.5) / total_variants`

#### 6. Curated SCORE_CONFIG (44 Traits)

Models were selected by maximizing chip coverage for reliability, with imputed coverage as secondary criteria:

| Category | Traits | Example PGS | Chip Coverage |
|----------|--------|-------------|---------------|
| Cardiovascular | CAD, Atrial Fib, Stroke, LDL, HDL, Triglycerides, PE | PGS000058 | 35–97% |
| Cancer | Breast, Colorectal, Prostate, Pancreatic, Testicular | PGS000028 | 51–97% |
| Metabolic | T2 Diabetes, BMI, WHR, HbA1c | PGS000125 | 62–100% |
| Neurological | Alzheimer's, Alzheimer's (APOE), Depression, ADHD | PGS000053 | 48–91% |
| Anthropometric | Height, Bone Density, Grip Strength, Lean Mass | PGS000297 | 11–93% |
| Autoimmune | Celiac, Crohn's, UC, Asthma | PGS000316 | 16–84% |
| Appearance | Skin Pigmentation, Hair Colour | PGS001244 | 14–17% |
| Behavioral | Risk-Taking, Cannabis Use, Nicotine Dependence | PGS001049 | 14–24% |
| Biomarkers | Uric Acid, Testosterone, CRP, IGF-1, Lp(a) | PGS000321 | 13–94% |

### Dual Coverage Reference

Every PRS result includes coverage metadata for transparency:

```json
{
  "coverage": {
    "matched": 200,
    "total": 204,
    "pct": 0.98,
    "chip_expected": 35.3,
    "imputed_raw": 98.0,
    "imputed_effective": 62.4
  }
}
```

### Ancestry-Aware Scoring

**Ancestry detection:** 14 ancestry-informative markers (AIMs) are used to classify the user's genetic ancestry as EUR, AFR, EAS, SAS, or AMR using Hardy-Weinberg log-likelihood scoring. Key markers include SLC24A5 (skin pigmentation), EDAR (hair thickness, EAS), and Duffy null (malaria resistance, AFR).

**Population-specific scoring:**
- Detected ancestry selects the appropriate reference population for both allele frequencies and distribution comparisons.
- Multi-ancestry allele frequencies are pre-computed from 1000 Genomes Phase 3: EUR (503 samples), AFR (661), EAS (504), SAS (489), AMR (347).
- For mixed ancestry (low confidence calls), percentiles are blended proportionally using softmax of AIM log-likelihoods across all populations.

**Ancestry-aware ML correction:**
- Ridge regression models include ancestry proportions as features: `[score, score², coverage, EUR%, AFR%, EAS%, SAS%, AMR%]`
- Trained on all 2,504 1000G samples (not EUR-only), so the model learns population-specific score adjustments.
- Per-population R² is tracked for each PGS model to assess cross-ancestry portability.

### ML Correction Pipeline (v3)

**Strand-ambiguous fix (v5 distributions):**
- Previous distributions (v4) skipped all A/T and C/G variants, but server scoring included them with direct matching.
- This caused systematic inflation for PGS with many strand-ambiguous variants (e.g., Height: 477/3290 = 14.5% ambiguous, 99.8% positive weights → z-score of 14+).
- v5 distributions now use the EXACT same scoring logic as the server, including direct matching for strand-ambiguous SNPs.

**Training pipeline:**
1. Score all 2,504 1000G samples across 22 chromosomes using harmonized PGS files with position-based matching.
2. Simulate chip-only scores by restricting to ~609K consumer chip positions.
3. Simulate imputed scores by adding ~20% noise to non-chip variant contributions.
4. Train per-PGS ancestry-aware Ridge regression models (44 models × chip + imputed).
5. Features: `[raw_score, raw_score², coverage, EUR%, AFR%, EAS%, SAS%, AMR%]`.
6. Cross-validated R² > 0.93 for imputed models; chip models variable (0.17–0.99) depending on coverage.

**Deployment:**
- Correction coefficients stored in `data/correction-models.json` (66KB).
- Applied at scoring time in Node.js — zero Python dependency at runtime.
- Backward-compatible: detects `ridge_ancestry` model type, falls back to legacy 3-coefficient model.

### Known Limitations

- PRS captures common variant risk only — rare pathogenic mutations are handled separately via ClinVar.
- Education/behavioral traits have low individual predictive power (population R² ~ 0.1).
- Hair color PGS can score high for blonde individuals on a "red hair" model (shared MC1R variants).
- PGS were predominantly developed in European cohorts; cross-ancestry performance varies by trait.
- Chip-only coverage limits accuracy for large PGS models (>5000 variants) — imputation strongly recommended.

---

## Pharmacogenomics (PGx)

### Dynamic Gene Coverage

34 pharmacogenes loaded dynamically from CPIC database (level A, A/B, and B evidence):

```
CYP2C9, CYP2C19, CYP2D6, CYP3A5, CYP4F2, DPYD, G6PD, NUDT15,
SLCO1B1, TPMT, UGT1A1, VKORC1, CYP2B6, CYP3A4, HLA-A, HLA-B, ...
```

### Star Allele Calling

CPIC allele definitions are loaded from the unified SQLite database (1,324 allele definitions across all pharmacogenes). For each gene:

1. **Load defining positions** from `cpic_sequence_locations` (rsID-mapped)
2. **Check user genotypes** at each defining position
3. **Score non-reference alleles** by matching defining variant patterns
4. **Determine zygosity** (homozygous vs heterozygous) from allele counts
5. **Fall back to reference** (*1) when no variant alleles match
6. **Look up diplotype phenotype** from `cpic_diplotypes` (110,346 entries)

Example output:
```
CYP2C9:  *1/*2  → Intermediate Metabolizer (6/6 positions covered)
CYP3A5:  *3/*3  → Poor Metabolizer (1/1 positions covered)
G6PD:    B/B    → Normal (3/3 positions covered)
```

### Clinical Data Integration

Three layers of pharmacogenomic data are provided per gene:

- **CPIC Recommendations:** Phenotype-specific dosing guidance with classification (Strong/Moderate)
- **CPIC Pairs:** Gene-drug interaction evidence levels
- **PharmGKB Annotations:** Variant-level drug-phenotype associations with significance scores

### Drug API

- `GET /api/drugs?q=warfarin` — Autocomplete search by drug name or gene
- `GET /api/drugs/:name` — Full detail with recommendations, gene pairs, and annotations

---

## Features

### DNA File Upload & Parsing
- **Supported formats:** 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, Nebula Genomics, VCF, zip archives
- Auto-format detection from file headers with filename-based fallback
- 100MB max file size, rate-limited (20 uploads/hour/IP)
- Beta access control via single-use codes (HELIX-XXXX-XXXX)

### Pathogenic Variant Scanning
- Full ClinVar database matching (400,000+ entries)
- Filters: pathogenic, likely pathogenic, risk factors, drug response
- SNPedia annotation matching (16,000+ entries) with gene, summary, and wiki content
- Deduplication across sources (ClinVar prioritized)

### Convergence Detection
- **9 disease pathway groups:** gastrointestinal cancer, digestive, cardiovascular, cardiometabolic, respiratory, neurodegenerative, mental health, autoimmune, kidney
- Intra-PRS convergence: flags when multiple scores in same pathway are elevated
- ClinVar+PRS convergence: monogenic pathogenic variants corroborated by polygenic risk
- Gene-to-category mapping (BRCA1/2, TP53, SCN5A, APOE, etc.)

### Deep Imputation (Optional)
- Beagle 5.5 imputation expands ~600K chip variants to ~28M
- 1000 Genomes Phase 3 reference panel, per-chromosome processing
- Imputed genotypes stored in SQLite with dosage and quality scores (gp_max)
- Quality tiers: High (GP >= 0.9), Medium (GP 0.7–0.9), Low (GP < 0.7)
- Typically ~50% high quality, ~25% medium, ~25% low for consumer chip data
- Progress tracking per chromosome with real-time status updates

### AI-Powered Protocols
- Claude generates personalized longevity protocols
- Input: top risks, strengths, convergence signals, ClinVar findings, PGx phenotypes, filtered variants
- Output: executive summary, domain summaries (cardiac, cancer, metabolism, brain, immune, drug, longevity), supplement recommendations (tiered), training protocols, lifespan interventions
- SNP-level annotations with category, risk, importance, and descriptions

### Reports & Export
- **Interactive HTML report** with Dashboard, Explorer, and Protocols tabs
- **4 PDF formats:**
  - Health Summary — essential PRS, risk stratification, key findings (~10-15 pages)
  - Doctor's Report — pharmacogenomics, drug-gene interactions, clinical format (~15-25 pages)
  - Full Analysis — all traits, complete ClinVar/SNPedia, convergence, protocols (~30-50 pages)
  - Health Protocol — actionable supplement and lifestyle interventions
- **JSON export** — machine-readable genetic profile for AI assistants
- **CSV export** — tabular data
- **Email delivery** via Resend API with all PDF attachments

### Variant Explorer
- Browse all matched ClinVar and SNPedia variants
- Side panel with gene name, genotype, significance, phenotype, related traits
- Searchable and filterable by source

### Trait Classification
- **Importance tiers:** essential (cancer, cardiovascular, metabolic, neurological, etc.), useful (biomarker, musculoskeletal, cognitive, etc.), technical
- **Polarity detection:** high-is-bad vs high-is-good per trait
- Auto-generated trait dictionary with friendly names, categories, and domain summaries

## Privacy & Security

### Zero Data Retention

- **Genetic data exists in volatile memory only** — deleted immediately after analysis completes
- Uploaded DNA files are unlinked from disk the moment scoring finishes (before the HTTP response)
- Reports, PDFs, and all intermediate files are automatically purged after 2 hours
- No user accounts, no cookies, no tracking, no IP logging
- Email used only for report delivery — not stored after send

### Cryptographic Deletion Certificates

Every user receives an automated **Data Deletion Certificate** when their data is purged:

- **SHA-256 attestation hash** — deterministic proof that the certificate was generated by the server at a specific time: `SHA256(certificateId:jobId:deletedAt:salt)`
- **PDF certificate** attached to the email — one-page document listing every data category destroyed
- **Append-only audit log** (`deletion-log.jsonl`) — tamper-evident record of all deletions (contains no genetic data, only masked emails and certificate hashes)

**Deleted items attested per certificate:**
1. Uploaded DNA genotype file
2. Imputation intermediate files (if applicable)
3. Genomic report (HTML, JSON data)
4. PRS scoring results
5. ClinVar & SNPedia match data
6. AI-generated health protocols
7. PDF reports
8. Job metadata and session data

**Deletion timeline:**
- DNA file: deleted immediately after analysis (~2 minutes)
- All other data: deleted 2 hours after upload
- Certificate email: sent automatically at deletion time
- Server restart recovery: stale jobs are cleaned up and certificates sent on boot

### Infrastructure Security

- Path traversal protection via UUID validation on all endpoints
- AI processing via local Claude CLI (no external data transmission)
- Rate limiting on uploads (20/hour/IP)
- File size limits (100MB max)
- Single concurrent analysis to prevent resource exhaustion

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| POST | `/api/validate-beta` | Validate a beta access code |
| POST | `/api/upload` | Upload DNA file (multipart form) |
| POST | `/api/analyze` | Start analysis (standard or imputed) |
| GET | `/api/status/:jobId` | Real-time job progress |
| GET | `/api/imputation-available` | Check if imputation is configured |
| GET | `/api/drugs?q=` | Drug autocomplete search |
| GET | `/api/drugs/:name` | Drug detail with PGx recommendations |
| GET | `/api/pdf/compact/:report` | Download health summary PDF |
| GET | `/api/pdf/doctor/:report` | Download doctor's report PDF |
| GET | `/api/pdf/full/:report` | Download full analysis PDF |
| GET | `/api/pdf/protocol/:report` | Download health protocol PDF |
| GET | `/reports/:jobId/index.html` | Interactive web report |

## Tech Stack

- **Runtime:** Node.js 22 + Express 4.21
- **Database:** SQLite3 (better-sqlite3) — unified DB with ClinVar, CPIC, PharmGKB, PGS metadata
- **PDF generation:** PDFKit
- **AI:** Anthropic Claude SDK
- **Imputation:** Beagle 5.5 (Java)
- **PRS reference:** 1000 Genomes Phase 3 (2,504 samples, 6 ancestry groups)
- **PGx reference:** CPIC (Clinical Pharmacogenetics Implementation Consortium)
- **Frontend:** Vanilla HTML/CSS/JS
- **ML correction:** Ancestry-aware Ridge regression (44 PGS models, 8-feature, trained on 2504 1000G samples)

## Data Architecture

```
data/
  helix-unified.db              # SQLite — ClinVar, SNPedia, CPIC tables, PGS metadata (~1.2GB)
  pgs-compact/                  # PGS scoring files: rsid\tea\toa\tweight (~450MB, 1200+ models)
  pgs-coverage-dual.json        # Chip + imputed coverage per PGS model
  pgs-chip-coverage.json        # Chip-only coverage reference
  correction-models.json        # Ancestry-aware Ridge regression coefficients (44 PGS × 8 features)
  prs-eur-freqs-1kg.json.gz     # EUR allele frequencies from 1000G (~27MB, 3.5M variants)
  prs-afr-freqs-1kg.json.gz     # AFR allele frequencies from 1000G (~27MB)
  prs-eas-freqs-1kg.json.gz     # EAS allele frequencies from 1000G (~26MB)
  prs-sas-freqs-1kg.json.gz     # SAS allele frequencies from 1000G (~27MB)
  prs-amr-freqs-1kg.json.gz     # AMR allele frequencies from 1000G (~26MB)
  prs-all-freqs-1kg.json.gz     # All-population allele frequencies (~29MB)

prs-distributions-v5.json.gz    # Population score distributions (2504 samples, 5 pops + ALL, 1202 PGS) (~39MB)

1000g/                          # 1000 Genomes Phase 3 VCFs (chr1-22) + panel.txt (~14GB)
pgs-harmonized/                 # Harmonized PGS files with chromosome positions (~927MB)
```

### SQLite Tables (helix-unified.db)

| Table | Rows | Description |
|-------|------|-------------|
| clinvar | 400K+ | ClinVar pathogenic/VUS variants |
| snpedia | 16K+ | SNPedia annotations |
| cpic_recommendations | 2,159 | CPIC dosing recommendations |
| cpic_pairs | 634 | Gene-drug interaction pairs |
| cpic_diplotypes | 110,346 | Diplotype-to-phenotype mappings |
| cpic_allele_definitions | 1,324 | Star allele definitions |
| cpic_sequence_locations | 1,184 | Allele-defining variant positions |
| cpic_allele_location_values | 3,091 | Variant alleles per allele definition |
| pharmgkb_annotations | 3,804 | PharmGKB variant-drug annotations |
| pgs_metadata | 1,200+ | PGS trait names and categories |

## Project Structure

```
server.js              # Express server & full analysis pipeline (~4500 lines)
package.json

public/                # Frontend (HTML, CSS, assets)
scripts/               # Runtime analysis & calibration scripts
tools/                 # Offline DB build, data processing, patching
infra/                 # Deployment (setup.sh) and backup (backup.sh)
docs/                  # VERSION.md, changelogs

# Gitignored (large data, not tracked):
data/                  # Reference databases (~37GB)
1000g/                 # 1000 Genomes Phase 3 reference panel (~14GB)
uploads/               # Temporary user DNA files
reports/               # Generated reports (auto-cleaned)
pgs-harmonized/        # Harmonized PGS scoring files
node_modules/
venv/
```

## Setup

```bash
# Clone and install
git clone <repo-url> /opt/helix-new
cd /opt/helix-new
npm install

# Configure
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY

# Reference data (not in git — obtain separately)
# Place helix-unified.db and PGS files in data/
# Place 1000 Genomes reference in 1000g/ (for imputation)

# Run
node --max-old-space-size=16384 server.js

# Or via systemd
sudo systemctl start helix
```

## License

Proprietary. All rights reserved.

