#!/usr/bin/env python3
"""Validate PGP top 50 genotype files: check chipset completeness + phenotype richness.
Flags files with <500K SNPs as incomplete."""

import os, json, subprocess, gzip, bz2

DATA_DIR = os.environ.get("HELIX_DATA_DIR", "./data")
GENO_DIR = os.path.join(DATA_DIR, "pgp-top50/genotypes/")
PHENO_FILE = os.path.join(DATA_DIR, "pgp-top50/pgp_top50_phenotypes.json")

with open(PHENO_FILE) as f:
    pheno = json.load(f)

results = []

for fn in sorted(os.listdir(GENO_DIR)):
    fpath = os.path.join(GENO_DIR, fn)

    # Extract hu_id
    import re
    hu_match = re.search(r'hu[A-F0-9]{6}', fn, re.IGNORECASE)
    hu_id = hu_match.group() if hu_match else "unknown"

    # Determine file type
    if "complete_genomics" in fn:
        chip_type = "Complete Genomics"
    elif "family_tree" in fn or "ftdna" in fn.lower():
        chip_type = "FTDNA"
    elif "23andme" in fn:
        chip_type = "23andMe"
    elif "ancestry" in fn.lower():
        chip_type = "AncestryDNA"
    else:
        chip_type = "Unknown"

    # Count lines (SNPs)
    line_count = 0
    try:
        if fn.endswith('.bz2'):
            with bz2.open(fpath, 'rt') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        line_count += 1
        elif fn.endswith('.gz'):
            with gzip.open(fpath, 'rt') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        line_count += 1
        elif fn.endswith('.zip'):
            import zipfile
            with zipfile.ZipFile(fpath) as z:
                for zf in z.namelist():
                    if zf.startswith('__') or zf.endswith('/'):
                        continue
                    with z.open(zf) as f:
                        for line in f:
                            l = line.decode('utf-8', errors='ignore')
                            if l.strip() and not l.startswith('#'):
                                line_count += 1
        else:
            with open(fpath) as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        line_count += 1
    except Exception as e:
        line_count = -1
        print(f"ERROR reading {fn}: {e}")

    # Get file size
    fsize = os.path.getsize(fpath)

    # Get phenotype data
    p = pheno.get(hu_id, {})
    profile = p.get("profile_data", {})
    traits = profile.get("traits", {})
    conditions = profile.get("conditions", [])
    n_surveys = p.get("survey_count", 0)
    n_responses = p.get("total_responses", 0)

    # Count medical-relevant traits
    medical_keywords = ["diabetes", "cancer", "heart", "cholesterol", "blood pressure",
                       "hypertension", "asthma", "depression", "anxiety", "arthritis",
                       "thyroid", "kidney", "liver", "weight", "height", "bmi",
                       "diagnosed", "medication", "allergy", "surgery"]
    medical_count = 0
    if isinstance(traits, dict):
        for k, v in traits.items():
            kl = (k + " " + str(v)).lower()
            if any(w in kl for w in medical_keywords):
                medical_count += 1

    # Has lab values?
    lab_keywords = ["cholesterol", "ldl", "hdl", "triglyceride", "glucose", "hba1c",
                   "creatinine", "bun", "ast", "alt", "tsh", "psa", "albumin"]
    has_labs = False
    if isinstance(traits, dict):
        for k in traits:
            if any(w in k.lower() for w in lab_keywords):
                has_labs = True
                break

    # Determine if chipset is full
    is_full = line_count >= 500000

    # Quality score: combines chip completeness + phenotype richness
    quality = 0
    if is_full:
        quality += 50
    elif line_count > 100000:
        quality += 25
    quality += min(n_responses / 20, 25)  # up to 25 for responses
    quality += min(medical_count / 5, 15)  # up to 15 for medical data
    if has_labs:
        quality += 10

    results.append({
        "rank": fn.split("_")[0],
        "hu_id": hu_id,
        "filename": fn,
        "chip_type": chip_type,
        "snp_count": line_count,
        "file_size_mb": round(fsize / 1024 / 1024, 1),
        "is_full_chip": is_full,
        "surveys": n_surveys,
        "responses": n_responses,
        "medical_traits": medical_count,
        "has_labs": has_labs,
        "conditions": len(conditions) if isinstance(conditions, list) else 0,
        "quality_score": round(quality, 1),
    })

# Sort by quality score
results.sort(key=lambda x: -x["quality_score"])

# Print summary
print("=" * 100)
print(f"{'Rank':<5} {'PGP ID':<10} {'Chip':<18} {'SNPs':>10} {'Full?':<6} {'Surveys':>8} {'Medical':>8} {'Labs':<5} {'Quality':>8}")
print("=" * 100)

full_count = 0
for r in results:
    status = "YES" if r["is_full_chip"] else "NO"
    labs = "YES" if r["has_labs"] else ""
    if r["is_full_chip"]:
        full_count += 1
    print(f"{r['rank']:<5} {r['hu_id']:<10} {r['chip_type']:<18} {r['snp_count']:>10,} {status:<6} {r['surveys']:>8} {r['medical_traits']:>8} {labs:<5} {r['quality_score']:>8.1f}")

print("=" * 100)
print(f"\nFull chipsets (>=500K SNPs): {full_count} / {len(results)}")
print(f"Ready for pipeline: {sum(1 for r in results if r['is_full_chip'] and r['quality_score'] >= 50)}")

# Save results
with open(os.path.join(DATA_DIR, "pgp-top50/chip_validation.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {os.path.join(DATA_DIR, pgp-top50/chip_validation.json)}")
