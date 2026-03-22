<p align="center">
  <img src="https://img.shields.io/badge/PGS_Models-3%2C550+-8b5cf6?style=for-the-badge" alt="PGS Models">
  <img src="https://img.shields.io/badge/ClinVar-400K+-06b6d4?style=for-the-badge" alt="ClinVar">
  <img src="https://img.shields.io/badge/Pharmacogenes-34_CPIC-10b981?style=for-the-badge" alt="CPIC">
  <img src="https://img.shields.io/badge/License-MIT_/_CC_BY_4.0-f59e0b?style=for-the-badge" alt="License">
</p>

<h1 align="center">Helix Open Research</h1>

<p align="center">
  <strong>Open research, open tools, open journey.</strong><br>
  <em>Genomics is too important to build behind closed doors.</em>
</p>

<p align="center">
  <a href="https://helixsequencing.com">Platform</a> ·
  <a href="journal/2026-03-prs-scoring-journey.md">Engineering Journal</a> ·
  <a href="tools/README.md">Research Tools</a> ·
  <a href="#get-involved">Collaborate</a>
</p>

---

## The Mission

We're building [Helix Sequencing](https://helixsequencing.com) — a privacy-first DNA analysis platform that turns a $99 consumer genotype file into 2,800+ polygenic risk scores, pharmacogenomics, pathogenic variant scanning, and AI-powered health protocols. Zero data retention.

Along the way, we've learned hard lessons. Spent money on experiments that didn't work. Hit walls nobody warned us about. **So we're sharing all of it** — the wins, the failures, the code, and the data.

> **Why?** Academic papers describe PRS methods in theory. Commercial platforms keep implementations proprietary. Nobody publishes the engineering reality — what actually works when you score 3,550 models against real consumer DNA data, what breaks, and what it costs. We believe this needs to change.

---

## Current Work: Iterative PRS Calibration (March 2026)

We are systematically processing **1,000+ genomes** from the [Personal Genome Project](https://www.personalgenomes.org/) to validate and improve polygenic risk score accuracy. PGP participants have consented to public sharing of both genetic data and medical records, making it the ideal ground-truth dataset.

### Methodology

Rather than scoring everything at once and hoping for the best, we're running an **iterative batch validation protocol**:

1. **Score** a batch of ~50 PGP participants through the full analysis pipeline
2. **Validate** PRS predictions against participants' documented medical diagnoses
3. **Identify** discordant cases — high-risk scores with no diagnosis, low-risk scores in diagnosed individuals
4. **Trace** discordances to specific PGS models and root-cause the failure mode
5. **Recalibrate** — remove underperforming models, verify ensemble composition
6. **Repeat** with the next batch, carrying forward all improvements

This prioritises stability over speed. Each batch builds on validated improvements from the previous cycle.

### What We've Found So Far

Validation against 64 PGP participants with medical records has identified several categories of problematic PGS models:

| Failure Mode | Example | Impact |
|:-------------|:--------|:-------|
| **Distribution bias** | Childhood-onset asthma model with population mean at 68th percentile instead of ~50th | Systematic false positives — nearly everyone scores elevated |
| **Wrong phenotype target** | Rheumatoid *factor* model scored as rheumatoid *arthritis* | Introduces noise from a related but distinct biomarker |
| **Duplicate models** | Two models from the same 2009 study, both included in ensemble | Overweights a single data source, reduces diversity |
| **Below polygenic threshold** | Models with <20 variants passing through individual scoring | Not truly polygenic — single-gene effects masquerading as risk scores |

Each batch cycle has improved concordance between predicted risk and actual diagnoses. The blocklist has grown from 0 to 54 models across multiple validation rounds.

### Roadmap

- **Ongoing:** Continue processing PGP batches (~900 remaining), validate after each cycle
- **Planned:** Add optional self-reported diagnosis collection to the platform, enabling continuous validation from every user who opts in
- **Planned:** Per-model discrimination metrics (AUROC against our own validated cohort) once sample sizes reach statistical significance (~200+ per trait)
- **Planned:** Ancestry-specific calibration — current validation is predominantly EUR; expanding to multi-ancestry cohorts

### Validation Metrics

| Metric | Current (n=64) | Target (n=500+) |
|:-------|:---------------|:-----------------|
| PRS-phenotype concordance | 51.4% | >70% |
| Blocklisted PGS models | 54 | Converging |
| Active ensemble traits | 479 | Stable |
| PGP genomes processed | 127 | 1,000+ |

---

## Engineering Journal

Honest post-mortems with real numbers, real costs, and real failures. Not polished marketing.

<table>
<tr>
<td width="100"><strong>Entry 1</strong><br><sub>Mar 2026</sub></td>
<td>
  <strong><a href="journal/2026-03-prs-scoring-journey.md">How We Built Accurate PRS from a $99 DNA Kit</a></strong><br>
  <sub>Scored 2,800+ disease risk models. Naive version was wildly wrong. Three weeks of "improvements" made it worse. Going back to basics produced the best results.</sub><br><br>
  <img src="https://img.shields.io/badge/Validation-4%2C257_samples-06b6d4?style=flat-square" alt="">
  <img src="https://img.shields.io/badge/Height_r-0.107_(p%3D0.044)-8b5cf6?style=flat-square" alt="">
  <img src="https://img.shields.io/badge/Build_Cost-%24400--800-10b981?style=flat-square" alt="">
</td>
</tr>
</table>

---

## Research Tools

Standalone Python scripts for PRS research. All work with publicly available data. **[Full documentation →](tools/README.md)**

| Category | Tools | What They Do |
|:---------|:------|:-------------|
| **Data Collection** | `scrape-pgp.py` `download-pgs.py` `fetch-ancestry.py` | Scrape PGP phenotypes, bulk download PGS Catalog models, fetch ancestry metadata |
| **Scoring** | `gpu-scorer.py` `build-distributions.py` | GPU-accelerated PRS via PyTorch sparse matrix multiplication, 1000 Genomes population distributions |
| **Validation** | `validate-opensnp.py` `validate-pgp.py` | Test against 4,257 OpenSNP samples and 943 PGP genomes with medical records |
| **QC** | `allele-alignment.py` `condition-mapper.py` | Strand-ambiguous SNP handling, clinical text-to-trait mapping (170+ patterns) |

---

## Platform at a Glance

<table>
<tr>
<td align="center" width="16%"><strong>3,550+</strong><br><sub>PGS Models</sub></td>
<td align="center" width="16%"><strong>400K+</strong><br><sub>ClinVar Variants</sub></td>
<td align="center" width="16%"><strong>34</strong><br><sub>CPIC Genes</sub></td>
<td align="center" width="16%"><strong>28M</strong><br><sub>Imputed Variants</sub></td>
<td align="center" width="16%"><strong>16K+</strong><br><sub>SNPedia Entries</sub></td>
<td align="center" width="16%"><strong>7</strong><br><sub>AI Agents</sub></td>
</tr>
</table>

**16 integrated databases** — PGS Catalog, ClinVar, MyVariant.info (CADD + gnomAD), MyDisease.info, MyChem.info, CPIC, DisGeNET, GWAS Catalog, Orphanet, PharmGKB, AlphaMissense, ClinGen, SNPedia, 1000 Genomes — unified in a 2GB SQLite database indexed by rsID.

**Full database documentation →** [HelixDocuments](https://github.com/HelixGenomics/HelixDocuments)

---

## Get Involved

We're actively looking for collaborators:

| Area | What We Need |
|:-----|:-------------|
| **PRS Validation** | Large-scale clinical dataset validation, real diagnosis ground truth |
| **Ancestry Accuracy** | Non-EUR PRS calibration, cross-ancestry portability research |
| **Rare Variants** | Combined rare + common variant scoring (RICE integration) |
| **Family Studies** | Trio analysis for rare disease — directly relevant to our Trisomy 9 research |
| **Pharmacogenomics** | Star allele calling improvements, CYP variant edge cases |

**Open an issue.** Submit a PR. Fork the tools and run your own validation. Or just use what's here and let us know what you find.

> *The more people working on this openly, the better the science gets for everyone.*

---

## License

| Content | License |
|:--------|:--------|
| Journal & documentation | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — share and adapt with attribution |
| Research tools | [MIT](https://opensource.org/licenses/MIT) — use freely for any purpose |

---

<p align="center">
  <a href="https://helixsequencing.com"><strong>helixsequencing.com</strong></a><br>
  <sub>Privacy-first DNA analysis · Zero data retention · Cryptographic deletion certificates</sub>
</p>
