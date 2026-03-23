#!/usr/bin/env python3
"""
Build ClinVar Star Rating Enrichment
=====================================

Parses the ClinVar VCF to extract review status star ratings for each variant,
then enriches an existing ClinVar database with star-weighted evidence scores.

ClinVar review stars (0-4) indicate the strength of evidence behind a classification:
  4 = Practice guideline
  3 = Expert panel reviewed
  2 = Multiple submitters, no conflict
  1 = Single submitter
  0 = No assertion criteria provided

Usage:
  python build-clinvar-stars.py --vcf clinvar.vcf.gz --db variants.db

The VCF is available from:
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

Output:
  Adds/updates a `clinvar_stars` table in the SQLite database with:
  - rsid, clinvar_id, review_status, star_count, significance, last_evaluated
"""

import argparse
import gzip
import re
import sqlite3
import sys
from datetime import datetime


STAR_MAP = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_for_the_single_variant": 0,
}

EVIDENCE_WEIGHTS = {4: 1.0, 3: 0.9, 2: 0.7, 1: 0.4, 0: 0.1}


def parse_info(info_str):
    """Parse VCF INFO field into dict."""
    fields = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            fields[key] = val
        else:
            fields[item] = True
    return fields


def parse_vcf(vcf_path):
    """Parse ClinVar VCF, yielding (rsid, clinvar_id, review_status, stars, significance)."""
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t", 8)
            if len(cols) < 8:
                continue

            chrom, pos, var_id, ref, alt, qual, filt, info_str = cols[:8]

            info = parse_info(info_str)

            # Extract RS number
            rsid = info.get("RS", "")
            if rsid:
                rsid = f"rs{rsid}"
            else:
                continue  # Skip variants without rsIDs

            clinvar_id = var_id
            review = info.get("CLNREVSTAT", "").lower().replace("_", "_")
            significance = info.get("CLNSIG", "").replace("_", " ")

            stars = 0
            for pattern, star_count in STAR_MAP.items():
                if pattern in review:
                    stars = star_count
                    break

            yield rsid, clinvar_id, review, stars, significance


def main():
    parser = argparse.ArgumentParser(description="Build ClinVar star rating enrichment")
    parser.add_argument("--vcf", required=True, help="Path to clinvar.vcf.gz")
    parser.add_argument("--db", required=True, help="Path to SQLite database")
    args = parser.parse_args()

    print(f"Parsing ClinVar VCF: {args.vcf}")
    print(f"Target database: {args.db}")

    conn = sqlite3.connect(args.db)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS clinvar_stars (
            rsid TEXT PRIMARY KEY,
            clinvar_id TEXT,
            review_status TEXT,
            star_count INTEGER,
            significance TEXT,
            evidence_weight REAL
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_clinvar_stars_rsid ON clinvar_stars(rsid)")

    batch = []
    count = 0
    for rsid, clinvar_id, review, stars, significance in parse_vcf(args.vcf):
        weight = EVIDENCE_WEIGHTS.get(stars, 0.1)
        batch.append((rsid, clinvar_id, review, stars, significance, weight))
        count += 1

        if len(batch) >= 10000:
            conn.executemany(
                "INSERT OR REPLACE INTO clinvar_stars VALUES (?, ?, ?, ?, ?, ?)",
                batch,
            )
            conn.commit()
            batch = []
            print(f"  Processed {count:,} variants...", end="\r")

    if batch:
        conn.executemany(
            "INSERT OR REPLACE INTO clinvar_stars VALUES (?, ?, ?, ?, ?, ?)",
            batch,
        )
        conn.commit()

    # Print summary
    cursor = conn.execute("SELECT star_count, COUNT(*) FROM clinvar_stars GROUP BY star_count ORDER BY star_count")
    print(f"\nDone. {count:,} variants processed.")
    print("\nStar distribution:")
    for stars, cnt in cursor:
        print(f"  {stars} stars: {cnt:,}")

    conn.close()


if __name__ == "__main__":
    main()
