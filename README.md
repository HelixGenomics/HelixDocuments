<p align="center">
  <img src="https://img.shields.io/badge/PGS_Models-3%2C550+-8b5cf6?style=for-the-badge" alt="PGS Models">
  <img src="https://img.shields.io/badge/ClinVar-400K+-06b6d4?style=for-the-badge" alt="ClinVar">
  <img src="https://img.shields.io/badge/Pharmacogenes-34_CPIC-10b981?style=for-the-badge" alt="CPIC">
  <img src="https://img.shields.io/badge/License-MIT_/_CC_BY_4.0-f59e0b?style=for-the-badge" alt="License">
</p>

<h1 align="center">🧬 Helix Open Research</h1>

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

## 📖 Engineering Journal

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

## 🛠️ Research Tools

Standalone Python scripts for PRS research. All work with publicly available data. **[Full documentation →](tools/README.md)**

| Category | Tools | What They Do |
|:---------|:------|:-------------|
| **Data Collection** | `scrape-pgp.py` `download-pgs.py` `fetch-ancestry.py` | Scrape PGP phenotypes, bulk download PGS Catalog models, fetch ancestry metadata |
| **Scoring** | `gpu-scorer.py` `build-distributions.py` | GPU-accelerated PRS via PyTorch sparse matrix multiplication, 1000 Genomes population distributions |
| **Validation** | `validate-opensnp.py` `validate-pgp.py` | Test against 4,257 OpenSNP samples and 943 PGP genomes with medical records |
| **QC** | `allele-alignment.py` `condition-mapper.py` | Strand-ambiguous SNP handling, clinical text-to-trait mapping (170+ patterns) |

---

## 📊 Platform at a Glance

<table>
<tr>
<td align="center" width="16%"><strong>3,550+</strong><br><sub>PGS Models</sub></td>
<td align="center" width="16%"><strong>400K+</strong><br><sub>ClinVar Variants</sub></td>
<td align="center" width="16%"><strong>34</strong><br><sub>CPIC Genes</sub></td>
<td align="center" width="16%"><strong>28M</strong><br><sub>Imputed Variants</sub></td>
<td align="center" width="16%"><strong>16K+</strong><br><sub>SNPedia Entries</sub></td>
<td align="center" width="16%"><strong>6</strong><br><sub>AI Agents</sub></td>
</tr>
</table>

**16 integrated databases** — PGS Catalog, ClinVar, MyVariant.info (CADD + gnomAD), MyDisease.info, MyChem.info, CPIC, DisGeNET, GWAS Catalog, Orphanet, PharmGKB, AlphaMissense, ClinGen, SNPedia, 1000 Genomes — unified in a 2GB SQLite database indexed by rsID.

**Full database documentation →** [HelixDocuments](https://github.com/HelixGenomics/HelixDocuments)

---

## 🤝 Get Involved

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

## 📄 License

| Content | License |
|:--------|:--------|
| Journal & documentation | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — share and adapt with attribution |
| Research tools | [MIT](https://opensource.org/licenses/MIT) — use freely for any purpose |

---

<p align="center">
  <a href="https://helixsequencing.com"><strong>helixsequencing.com</strong></a><br>
  <sub>Privacy-first DNA analysis · Zero data retention · Cryptographic deletion certificates</sub>
</p>
