#!/usr/bin/env python3
"""
Health Index Simulation & Calibration Tool
============================================

Simulates the composite health index score across a cohort to calibrate
weighting parameters. The health index combines PRS risk scores, pathogenic
variant counts, carrier status, and pharmacogenomic actionability into a
single 0-100 score.

Usage:
  python sim-health-index.py --reports-dir /path/to/reports/ [--output results.csv]

The health index formula:
  score = 100 - (prs_penalty + variant_penalty + carrier_penalty - pharma_bonus)

Where:
  prs_penalty:     Weighted sum of high-risk PRS percentiles (>75th)
  variant_penalty: Count of pathogenic/likely pathogenic variants * severity weight
  carrier_penalty: Count of carrier variants * 0.5 (informational, not disease-causing)
  pharma_bonus:    Actionable pharmacogenomic findings (positive — knowledge is power)

The goal is calibration: the median score across a healthy population should be
~65-75, with genuinely high-risk individuals scoring 30-50 and very low-risk
scoring 85-95.
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict


# Example weights — adjust these to calibrate for your cohort
# These are starting-point values, not production weights
WEIGHTS = {
    "prs_high_risk_threshold": 75,      # Percentile above which PRS contributes penalty
    "prs_very_high_threshold": 90,      # Percentile for "very high risk"
    "prs_penalty_per_high": 0.5,        # Penalty per high-risk trait (tune this)
    "prs_penalty_per_very_high": 1.0,   # Penalty per very-high-risk trait (tune this)
    "prs_max_penalty": 30,              # Cap on total PRS penalty
    "variant_pathogenic_weight": 2.0,   # Per pathogenic variant
    "variant_likely_path_weight": 1.0,  # Per likely pathogenic variant
    "variant_max_penalty": 20,          # Cap on variant penalty
    "carrier_weight": 0.5,             # Per carrier variant
    "carrier_max_penalty": 5,           # Cap on carrier penalty
    "pharma_actionable_bonus": 1.0,     # Per actionable PGx finding
    "pharma_max_bonus": 10,             # Cap on PGx bonus
}


def load_prs(report_dir):
    """Load PRS scores from a report directory."""
    for fname in ["ensembled-prs.json", "prs-results.json"]:
        path = os.path.join(report_dir, fname)
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            if isinstance(data, list):
                return data
            elif isinstance(data, dict) and "scores" in data:
                return data["scores"]
    return []


def load_variants(report_dir):
    """Load variant data from explorer-data.json."""
    path = os.path.join(report_dir, "explorer-data.json")
    if not os.path.exists(path):
        return {"pathogenic": 0, "likely_pathogenic": 0, "carriers": 0}

    with open(path) as f:
        data = json.load(f)

    counts = {"pathogenic": 0, "likely_pathogenic": 0, "carriers": 0}
    variants = data if isinstance(data, list) else data.get("variants", [])
    for v in variants:
        sig = v.get("significance", "").lower()
        if "pathogenic" in sig and "likely" not in sig:
            counts["pathogenic"] += 1
        elif "likely pathogenic" in sig or "likely_pathogenic" in sig:
            counts["likely_pathogenic"] += 1

        if v.get("carrier", False) or v.get("zygosity") == "heterozygous":
            counts["carriers"] += 1

    return counts


def load_pharma(report_dir):
    """Load pharmacogenomic actionability count."""
    for fname in ["agent-input-cpic.json", "agent-focus-pharma.json"]:
        path = os.path.join(report_dir, fname)
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            if isinstance(data, list):
                return sum(1 for item in data if item.get("actionable", False))
            elif isinstance(data, dict):
                items = data.get("results", data.get("phenotypes", []))
                if isinstance(items, list):
                    return len([i for i in items if i.get("actionable", True)])
    return 0


def calculate_health_index(prs_scores, variant_counts, pharma_count, weights=WEIGHTS):
    """Calculate composite health index score."""
    # PRS penalty
    high_risk = sum(1 for s in prs_scores
                    if s.get("percentile", 50) >= weights["prs_high_risk_threshold"]
                    and s.get("high_is_bad", True))
    very_high = sum(1 for s in prs_scores
                    if s.get("percentile", 50) >= weights["prs_very_high_threshold"]
                    and s.get("high_is_bad", True))

    prs_penalty = min(
        high_risk * weights["prs_penalty_per_high"] +
        very_high * weights["prs_penalty_per_very_high"],
        weights["prs_max_penalty"]
    )

    # Variant penalty
    variant_penalty = min(
        variant_counts["pathogenic"] * weights["variant_pathogenic_weight"] +
        variant_counts["likely_pathogenic"] * weights["variant_likely_path_weight"],
        weights["variant_max_penalty"]
    )

    # Carrier penalty (informational)
    carrier_penalty = min(
        variant_counts["carriers"] * weights["carrier_weight"],
        weights["carrier_max_penalty"]
    )

    # Pharma bonus
    pharma_bonus = min(
        pharma_count * weights["pharma_actionable_bonus"],
        weights["pharma_max_bonus"]
    )

    score = 100 - prs_penalty - variant_penalty - carrier_penalty + pharma_bonus
    score = max(0, min(100, score))

    return {
        "score": round(score, 1),
        "prs_penalty": round(prs_penalty, 1),
        "variant_penalty": round(variant_penalty, 1),
        "carrier_penalty": round(carrier_penalty, 1),
        "pharma_bonus": round(pharma_bonus, 1),
        "n_high_risk_prs": high_risk,
        "n_very_high_prs": very_high,
        "n_pathogenic": variant_counts["pathogenic"],
        "n_likely_pathogenic": variant_counts["likely_pathogenic"],
        "n_carriers": variant_counts["carriers"],
        "n_pharma_actionable": pharma_count,
    }


def main():
    parser = argparse.ArgumentParser(description="Simulate health index across a cohort")
    parser.add_argument("--reports-dir", required=True, help="Directory containing report subdirectories")
    parser.add_argument("--output", default="health-index-simulation.csv", help="Output CSV path")
    args = parser.parse_args()

    results = []
    report_dirs = sorted(os.listdir(args.reports_dir))
    print(f"Processing {len(report_dirs)} reports...")

    for i, report_id in enumerate(report_dirs):
        report_dir = os.path.join(args.reports_dir, report_id)
        if not os.path.isdir(report_dir):
            continue

        prs = load_prs(report_dir)
        variants = load_variants(report_dir)
        pharma = load_pharma(report_dir)

        if not prs:
            continue

        index = calculate_health_index(prs, variants, pharma)
        index["report_id"] = report_id
        results.append(index)

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(report_dirs)}...")

    if not results:
        print("No reports with PRS data found.")
        sys.exit(1)

    # Write CSV
    fieldnames = ["report_id", "score", "prs_penalty", "variant_penalty",
                  "carrier_penalty", "pharma_bonus", "n_high_risk_prs",
                  "n_very_high_prs", "n_pathogenic", "n_likely_pathogenic",
                  "n_carriers", "n_pharma_actionable"]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    # Print summary stats
    scores = [r["score"] for r in results]
    scores.sort()
    mean_score = sum(scores) / len(scores)
    median_score = scores[len(scores) // 2]

    print(f"\nResults written to: {args.output}")
    print(f"Reports processed: {len(results)}")
    print(f"\nHealth Index Distribution:")
    print(f"  Mean:   {mean_score:.1f}")
    print(f"  Median: {median_score:.1f}")
    print(f"  Min:    {min(scores):.1f}")
    print(f"  Max:    {max(scores):.1f}")
    print(f"  Std:    {(sum((s - mean_score) ** 2 for s in scores) / len(scores)) ** 0.5:.1f}")

    # Distribution buckets
    buckets = {"90-100": 0, "75-89": 0, "60-74": 0, "45-59": 0, "30-44": 0, "0-29": 0}
    for s in scores:
        if s >= 90: buckets["90-100"] += 1
        elif s >= 75: buckets["75-89"] += 1
        elif s >= 60: buckets["60-74"] += 1
        elif s >= 45: buckets["45-59"] += 1
        elif s >= 30: buckets["30-44"] += 1
        else: buckets["0-29"] += 1

    print(f"\n  Score Range  | Count | Pct")
    print(f"  -----------  | ----- | ---")
    for bucket, count in buckets.items():
        pct = count / len(scores) * 100
        print(f"  {bucket:>10s}   | {count:5d} | {pct:.0f}%")

    # Calibration check
    if median_score < 60:
        print(f"\n  WARNING: Median too low ({median_score:.0f}). Consider reducing penalty weights.")
    elif median_score > 80:
        print(f"\n  WARNING: Median too high ({median_score:.0f}). Consider increasing penalty weights.")
    else:
        print(f"\n  Calibration looks good (median in 60-80 range).")


if __name__ == "__main__":
    main()
