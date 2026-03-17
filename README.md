# Helix Open Research

**Open research, open tools, open journey.**

This repository exists because genomics is too important to build behind closed doors.

We're building [Helix Sequencing](https://helixsequencing.com) — a privacy-first DNA analysis platform that turns consumer genotype files into actionable health insights using 3,550+ polygenic risk scores, pharmacogenomics, pathogenic variant scanning, and AI-powered longevity protocols.

Along the way, we've learned hard lessons about PRS accuracy, spent money on experiments that didn't work, hit walls that nobody warned us about, and figured things out that we think other people should know. So we're sharing all of it — the wins, the failures, the code, and the data.

## Why Open?

Polygenic risk scoring sits at the intersection of genomics and longevity science. The potential is enormous: identify your genetic predispositions before they manifest, intervene early, live longer and healthier. But the field is fragmented. Academic papers describe methods in theory. Commercial platforms keep their implementations proprietary. Nobody publishes the engineering reality — what actually works when you try to score 3,550 models against real consumer DNA data, what breaks, and what it costs.

We believe this needs to change.

The promise of genomics-informed longevity won't be realized by any single team working in isolation. It requires collective effort — researchers, engineers, clinicians, and individuals all contributing to a shared understanding of what genetic risk scores can and cannot tell us. The faster we share what we learn, the faster everyone benefits.

This repository is our contribution to that effort. Everything here is published with the hope that it saves someone else time, money, or frustration — and that it moves the entire field forward, even if only by a small step.

## What's Here

### Engineering Journal

Transparent, honest documentation of our R&D process. Not polished marketing — real engineering post-mortems with actual numbers, actual costs, and actual failures.

| Entry | Date | Topic |
|-------|------|-------|
| [Building Accurate PRS from Consumer DNA Data](journal/2026-03-prs-scoring-journey.md) | March 2026 | Five phases of PRS pipeline development. Validation results (r=0.107 for height). The PRS-CSx experiment that didn't improve accuracy. A 95% PGP imputation failure rate. Why correction models pulled everything to the mean. Lessons learned and roadmap. |

### Research Tools

Standalone Python scripts for PRS research, all working with publicly available data. [Full tool documentation](tools/README.md).

**Data collection** — Scrape PGP phenotypes, bulk download PGS Catalog models, fetch ancestry metadata

**Scoring** — GPU-accelerated PRS scoring via PyTorch sparse matrix multiplication, 1000 Genomes population distribution building, full batch pipelines for Vast.ai

**Validation** — Validate PRS against OpenSNP phenotypes (4,257 samples) and PGP medical records (943 profiles), select best model per trait per ancestry

**QC** — Allele alignment checking, clinical text-to-trait mapping with 170+ condition patterns

## About Helix Sequencing

Privacy-first DNA analysis. Upload a consumer genotype file (23andMe, AncestryDNA, MyHeritage, or VCF), receive a comprehensive genetic health report. Zero data retention — all user data is automatically deleted after analysis, with a cryptographic deletion certificate.

| What We Score | Scale |
|---------------|-------|
| Polygenic risk scores | 3,550+ models across cancer, cardiovascular, metabolic, neurological, autoimmune, and more |
| Pathogenic variants | 400,000+ ClinVar entries scanned |
| Pharmacogenomics | 34 CPIC genes, 110K+ diplotype-to-phenotype mappings |
| Variant annotations | 16,000+ SNPedia entries |
| Imputed variants | ~28M via Beagle 5.5 |
| AI health protocols | 6 parallel domain-specialized agents |

## Get Involved

If you're working on PRS accuracy, ancestry-aware scoring, consumer genomics, or longevity science — we'd love to hear from you. Open an issue, submit a PR, or just use the tools and let us know what you find.

The more people working on this openly, the better the science gets for everyone.

## License

Journal content and documentation: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — share and adapt freely with attribution.

Research tools: MIT — use freely for any purpose.
