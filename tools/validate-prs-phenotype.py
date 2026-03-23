#!/usr/bin/env python3
"""
PRS Validation with Sensitivity, Specificity, and AUC
=====================================================

Validates PRS scores against self-reported phenotype data from PGP participants.
Compares PRS percentile distributions between participants WITH vs WITHOUT each condition.

Data sources:
  - Phenotypes: /opt/helix/data/pgp-{top50,batch3,batch4}/pgp_top50_phenotypes.json
  - PRS scores: /opt/helix/reports/{jobId}/ensembled-prs.json
  - Job mapping: /opt/helix/data/master-batch-v3.log (COMPLETE lines)

Output:
  - /opt/helix/data/validation/prs-validation-results.csv
  - /opt/helix/data/validation/prs-validation-summary.txt
"""

import json
import os
import re
import csv
import sys
from collections import defaultdict
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────────────────────

PHENOTYPE_FILES = [
    "/opt/helix/data/pgp-top50/pgp_top50_phenotypes.json",
    "/opt/helix/data/pgp-batch3/pgp_top50_phenotypes.json",
    "/opt/helix/data/pgp-batch4/pgp_top50_phenotypes.json",
]
REPORTS_DIR = "/opt/helix/reports"
BATCH_LOGS = [
    "/opt/helix/data/master-batch-v3.log",  # Only v3 — earlier batches had unfixed bugs
]
OUTPUT_DIR = "/opt/helix/data/validation"
MIN_GROUP_SIZE = 3  # Minimum participants per group for meaningful comparison

# ── Condition → PRS trait mapping ──────────────────────────────────────────────
# Maps self-reported conditions (from surveys/traits) to PRS trait names
# PRS traits are lowercase, normalized

CONDITION_MAP = {
    # Cardiovascular
    "Hypertension": "hypertension",
    "High cholesterol (hypercholesterolemia)": "cholesterol",
    "High triglycerides (hypertriglyceridemia)": "triglycerides",
    "Cardiac arrhythmia": "atrial fibrillation",
    "Atrial fibrillation": "atrial fibrillation",
    "Coronary artery disease": "coronary artery disease",
    "Heart attack": "coronary artery disease",
    "Mitral valve prolapse": "mitral valve prolapse",
    "Hemorrhoids": "hemorrhoids",
    "Deep vein thrombosis": "blood clot or deep vein thrombosis",
    "Varicose veins": "varicose veins",

    # Metabolic
    "Type 2 Diabetes": "type 2 diabetes",
    "Type 1 Diabetes": "type 1 diabetes",
    "Gout": "gout",
    "Hypothyroidism": "hypothyroidism",
    "Hyperthyroidism": "hyperthyroidism",
    "Iron deficiency anemia": "iron deficiency anemia",

    # Respiratory
    "Asthma": "asthma",
    "Asthma (Adult)": "asthma",
    "Asthma (Childhood)": "asthma",
    "Chronic obstructive pulmonary disease (COPD": "chronic obstructive pulmonary disease",
    "Chronic bronchitis": "chronic obstructive pulmonary disease",
    "Emphysema": "chronic obstructive pulmonary disease",
    "Pneumonia": "pneumonia",
    "Allergic rhinitis": "allergic rhinitis",

    # Musculoskeletal
    "Osteoarthritis": "osteoarthritis",
    "Osteoporosis": "osteoporosis",
    "Rheumatoid arthritis": "rheumatoid arthritis",
    "Scoliosis": "scoliosis",
    "Sciatica": "sciatica",
    "Plantar fasciitis": "plantar fasciitis",
    "Carpal tunnel syndrome": "carpal tunnel syndrome",
    "Dupuytren's contracture": "contracture of palmar fascia dupuytrens disease",
    "Bunions": "hallux valgus",
    "Tennis elbow": "tennis elbow",
    "Ankylosing spondylitis": "ankylosing spondylitis",

    # Gastrointestinal
    "Gastroesophageal reflux disease (GERD)": "gastroesophageal reflux disease",
    "Irritable bowel syndrome (IBS)": "irritable bowel syndrome",
    "Celiac disease": "celiac disease",
    "Diverticulosis": "diverticular disease",
    "Colon polyps": "benign neoplasm of colon",
    "Inguinal hernia": "inguinal hernia",
    "Hiatal hernia": "hiatal hernia",
    "Gallstones": "cholelithiasis",
    "Cholelithiasis": "cholelithiasis",

    # Dermatological
    "Eczema": "atopic eczema",
    "Acne": "acne",
    "Rosacea": "rosacea",
    "Psoriasis": "psoriasis",
    "Skin tags": "skin tags",
    "Non-melanoma skin cancer": "basal cell carcinoma",
    "Melanoma": "melanoma",
    "Dandruff": "dandruff",

    # Cancer
    "Breast cancer": "breast cancer",
    "Prostate cancer": "prostate cancer",
    "Colorectal cancer": "colorectal cancer",
    "Lung cancer": "lung cancer",
    "Thyroid cancer": "thyroid cancer",
    "Bladder cancer": "bladder cancer",

    # Neurological
    "Migraine without aura": "migraine",
    "Migraine with aura": "migraine",
    "Tinnitus": "tinnitus",
    "Age-related hearing loss": "hearing difficulty",

    # Eye
    "Myopia (Nearsightedness)": "myopia",
    "Age-related macular degeneration": "age related macular degeneration",
    "Age-related cataract": "cataract",
    "Glaucoma": "glaucoma",

    # Kidney/Urinary
    "Kidney stones": "calculus of kidney and ureter",
    "Chronic kidney disease": "chronic kidney disease",

    # Psychiatric
    "Depression": "major depressive disorder",
    "Anxiety": "generalised anxiety disorder",
    "ADHD": "attention deficit hyperactivity disorder",
    "Bipolar disorder": "bipolar disorder",
    "Schizophrenia": "schizophrenia",

    # Other
    "Dental cavities": "dental caries",
    "Chronic sinusitis": "chronic sinusitis",
    "Sleep apnea": "sleep apnoea",
    "Hypothyroidism": "hypothyroidism",
    "Thyroid nodule(s)": "thyroid nodule",
    "Chronic tonsillitis": "chronic tonsillitis",
    "Hair loss (includes female and male pattern baldness)": "male pattern baldness",
    "Endometriosis": "endometriosis",
    "Ovarian cysts": "ovarian cyst",
    "Uterine fibroids": "uterine fibroids",
    "Benign prostatic hypertrophy (BPH)": "benign prostatic hyperplasia",

    # Profile traits (date-format conditions)
    "Benign Prostatic Hypertrophy (BPH)": "benign prostatic hyperplasia",
    "Heart murmur": "heart murmur",
    "High Cholesterol": "cholesterol",
    "Thyroid Nodule": "thyroid nodule",
}


def load_phenotypes():
    """Load phenotypes from all batch files, return dict of hu_id -> set of conditions."""
    participants = {}

    for path in PHENOTYPE_FILES:
        if not os.path.exists(path):
            print(f"  Skipping {path} (not found)")
            continue

        with open(path) as f:
            data = json.load(f)

        if not data:
            continue

        print(f"  Loaded {path}: {len(data)} participants")

        for hu_id, entry in data.items():
            if hu_id not in participants:
                participants[hu_id] = {
                    "conditions": set(),
                    "no_conditions": set(),  # Explicitly denied conditions (from s23)
                }

            profile = entry.get("profile_data", {})
            survey = entry.get("survey_traits", {})

            # 1. Extract from profile_data.traits (date-format = diagnosed)
            traits = profile.get("traits", {})
            for trait_name, trait_val in traits.items():
                if trait_name in CONDITION_MAP:
                    participants[hu_id]["conditions"].add(CONDITION_MAP[trait_name])

            # 2. Extract from survey_traits - "Have you ever been diagnosed..." (comma-separated list)
            for key, val in survey.items():
                if "diagnosed" not in key.lower():
                    continue

                # Format 1: s23 bracket format [Condition Name]: Yes/No
                if "[" in key:
                    cond_name = key.split("[")[1].rstrip("]").strip()
                    val_str = str(val).strip().lower()
                    if cond_name in CONDITION_MAP:
                        prs_trait = CONDITION_MAP[cond_name]
                        if val_str == "yes":
                            participants[hu_id]["conditions"].add(prs_trait)
                        elif val_str == "no":
                            participants[hu_id]["no_conditions"].add(prs_trait)
                    continue

                # Format 2: Comma-separated list of conditions
                if val and str(val).strip().lower() not in ("no", "none", "n/a", ""):
                    for cond in str(val).split(","):
                        cond = cond.strip()
                        if cond in CONDITION_MAP:
                            participants[hu_id]["conditions"].add(CONDITION_MAP[cond])

    return participants


def load_job_mapping():
    """Parse all batch logs to get hu_id -> job_id mapping. Later logs override earlier ones."""
    mapping = {}
    for log_path in BATCH_LOGS:
        if not os.path.exists(log_path):
            print(f"  Skipping {log_path} (not found)")
            continue

        count = 0
        with open(log_path) as f:
            for line in f:
                m = re.search(r"COMPLETE \((\w+)\).*Job ([a-f0-9-]+)", line)
                if m:
                    mapping[m.group(1)] = m.group(2)
                    count += 1
        print(f"  {log_path}: {count} COMPLETE entries")

    print(f"  Total unique participants with jobs: {len(mapping)}")
    return mapping


def load_prs_scores(job_id):
    """Load PRS scores for a job. Supports both ensembled-prs.json (new) and prs-results.json (old)."""
    report_dir = os.path.join(REPORTS_DIR, job_id)

    # Try new format first
    path = os.path.join(report_dir, "ensembled-prs.json")
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
        scores = {}
        for item in data:
            scores[item["trait"]] = {
                "percentile": item["percentile"],
                "z_score": item.get("z_score", 0),
                "confidence": item.get("confidence", "unknown"),
                "high_is_bad": item.get("high_is_bad", True),
                "pgs_id": item.get("pgs_id", ""),
                "n_models": item.get("n_models", 1),
            }
        return scores

    # Fall back to old format
    path = os.path.join(report_dir, "prs-results.json")
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
        items = data.get("scores", data) if isinstance(data, dict) else data
        scores = {}
        for item in items:
            trait = item.get("trait", "")
            if not trait:
                continue
            scores[trait] = {
                "percentile": item.get("percentile", 50),
                "z_score": item.get("z_score", 0),
                "confidence": item.get("confidence", "unknown"),
                "high_is_bad": item.get("high_is_bad", True),
                "pgs_id": item.get("pgs_id", ""),
                "n_models": item.get("n_models", 1),
            }
        return scores

    return None


def calculate_auc(has_scores, no_scores):
    """Calculate AUC using Mann-Whitney U statistic (equivalent to AUC for binary classification)."""
    n_pos = len(has_scores)
    n_neg = len(no_scores)

    if n_pos == 0 or n_neg == 0:
        return None

    # Count concordant pairs: how often a positive case scores higher than a negative case
    concordant = 0
    tied = 0
    for p in has_scores:
        for n in no_scores:
            if p > n:
                concordant += 1
            elif p == n:
                tied += 0.5

    auc = (concordant + tied) / (n_pos * n_neg)
    return auc


def calculate_metrics(has_scores, no_scores, threshold=50.0):
    """Calculate sensitivity, specificity, and related metrics at a given percentile threshold."""
    # For high_is_bad traits: higher percentile = higher risk
    # Sensitivity: fraction of HAS_CONDITION scoring above threshold
    # Specificity: fraction of NO_CONDITION scoring below threshold

    tp = sum(1 for s in has_scores if s >= threshold)
    fn = sum(1 for s in has_scores if s < threshold)
    tn = sum(1 for s in no_scores if s < threshold)
    fp = sum(1 for s in no_scores if s >= threshold)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0

    return {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": ppv,
        "npv": npv,
        "tp": tp,
        "fn": fn,
        "tn": tn,
        "fp": fp,
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 70)
    print("PRS Validation with Specificity & AUC")
    print(f"Run: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Load data
    print("\n[1/4] Loading phenotypes...")
    phenotypes = load_phenotypes()
    print(f"  Total participants with phenotype data: {len(phenotypes)}")

    print("\n[2/4] Loading job mapping...")
    job_map = load_job_mapping()

    # Find participants with both phenotype data AND PRS scores
    print("\n[3/4] Loading PRS scores...")
    matched = {}
    for hu_id, pheno in phenotypes.items():
        if hu_id in job_map:
            job_id = job_map[hu_id]
            scores = load_prs_scores(job_id)
            if scores:
                matched[hu_id] = {"pheno": pheno, "scores": scores, "job_id": job_id}

    print(f"  Matched participants (phenotype + PRS): {len(matched)}")

    if not matched:
        print("ERROR: No matched participants found. Check job mapping and report paths.")
        sys.exit(1)

    # Build per-trait validation
    print("\n[4/4] Calculating per-condition metrics...")

    # Collect all PRS traits that have matching conditions
    prs_traits_with_data = set()
    for hu_id, data in matched.items():
        prs_traits_with_data.update(data["scores"].keys())

    condition_prs_traits = set(CONDITION_MAP.values())
    validatable = condition_prs_traits & prs_traits_with_data
    print(f"  PRS traits with condition mapping: {len(validatable)}")

    results = []

    for prs_trait in sorted(validatable):
        has_condition_scores = []
        no_condition_scores = []

        for hu_id, data in matched.items():
            if prs_trait not in data["scores"]:
                continue

            score_info = data["scores"][prs_trait]
            percentile = score_info["percentile"]
            high_is_bad = score_info.get("high_is_bad", True)

            has_it = prs_trait in data["pheno"]["conditions"]
            explicitly_no = prs_trait in data["pheno"]["no_conditions"]

            if has_it:
                has_condition_scores.append(percentile)
            elif explicitly_no or len(data["pheno"]["conditions"]) > 0:
                # Only count as "no condition" if they explicitly said no,
                # OR if they reported other conditions (meaning they engaged with the survey)
                no_condition_scores.append(percentile)

        n_has = len(has_condition_scores)
        n_no = len(no_condition_scores)

        if n_has < MIN_GROUP_SIZE or n_no < MIN_GROUP_SIZE:
            continue

        # Get score metadata from first available score
        sample_score = None
        for hu_id, data in matched.items():
            if prs_trait in data["scores"]:
                sample_score = data["scores"][prs_trait]
                break

        high_is_bad = sample_score.get("high_is_bad", True) if sample_score else True

        # If high_is_bad, higher percentile = higher risk = should be higher for affected
        # If NOT high_is_bad, we flip: affected should have LOWER percentile
        if high_is_bad:
            auc = calculate_auc(has_condition_scores, no_condition_scores)
            metrics = calculate_metrics(has_condition_scores, no_condition_scores, threshold=50.0)
        else:
            # Flip: affected should score lower
            auc = calculate_auc(
                [100 - s for s in has_condition_scores],
                [100 - s for s in no_condition_scores],
            )
            metrics = calculate_metrics(
                [100 - s for s in has_condition_scores],
                [100 - s for s in no_condition_scores],
                threshold=50.0,
            )

        mean_has = sum(has_condition_scores) / n_has
        mean_no = sum(no_condition_scores) / n_no
        separation = mean_has - mean_no if high_is_bad else mean_no - mean_has

        result = {
            "trait": prs_trait,
            "pgs_id": sample_score.get("pgs_id", "") if sample_score else "",
            "n_models": sample_score.get("n_models", 1) if sample_score else 1,
            "n_has": n_has,
            "n_no": n_no,
            "mean_has": round(mean_has, 1),
            "mean_no": round(mean_no, 1),
            "separation": round(separation, 1),
            "auc": round(auc, 3) if auc else None,
            "sensitivity": round(metrics["sensitivity"], 3),
            "specificity": round(metrics["specificity"], 3),
            "ppv": round(metrics["ppv"], 3),
            "npv": round(metrics["npv"], 3),
            "high_is_bad": high_is_bad,
            "confidence": sample_score.get("confidence", "") if sample_score else "",
            "flag": "",
        }

        # Flag problematic models
        if auc is not None:
            if auc < 0.50:
                result["flag"] = "WORSE_THAN_RANDOM"
            elif auc < 0.55:
                result["flag"] = "NO_BETTER_THAN_RANDOM"
            elif auc >= 0.65:
                result["flag"] = "GOOD"

        results.append(result)

    # Sort by AUC descending
    results.sort(key=lambda x: x["auc"] if x["auc"] is not None else 0, reverse=True)

    # Write CSV
    csv_path = os.path.join(OUTPUT_DIR, "prs-validation-results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "trait", "pgs_id", "n_models", "n_has", "n_no",
            "mean_has", "mean_no", "separation", "auc",
            "sensitivity", "specificity", "ppv", "npv",
            "high_is_bad", "confidence", "flag",
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n  Results written to: {csv_path}")

    # Write summary
    summary_path = os.path.join(OUTPUT_DIR, "prs-validation-summary.txt")
    with open(summary_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("PRS VALIDATION SUMMARY\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Participants matched: {len(matched)}\n")
        f.write(f"Conditions validated: {len(results)}\n")
        f.write("=" * 70 + "\n\n")

        # Summary stats
        aucs = [r["auc"] for r in results if r["auc"] is not None]
        if aucs:
            f.write(f"Mean AUC: {sum(aucs)/len(aucs):.3f}\n")
            f.write(f"Median AUC: {sorted(aucs)[len(aucs)//2]:.3f}\n")
            f.write(f"AUC ≥ 0.65 (good): {sum(1 for a in aucs if a >= 0.65)}\n")
            f.write(f"AUC 0.55-0.65 (moderate): {sum(1 for a in aucs if 0.55 <= a < 0.65)}\n")
            f.write(f"AUC < 0.55 (poor): {sum(1 for a in aucs if a < 0.55)}\n")
            f.write(f"AUC < 0.50 (worse than random): {sum(1 for a in aucs if a < 0.50)}\n")

        # Best performing
        f.write("\n── TOP PERFORMING (AUC ≥ 0.60) " + "─" * 40 + "\n")
        for r in results:
            if r["auc"] and r["auc"] >= 0.60:
                f.write(f"  {r['trait']:45s} AUC={r['auc']:.3f}  has={r['n_has']:3d}  no={r['n_no']:3d}  "
                        f"mean_diff={r['separation']:+5.1f}  sens={r['sensitivity']:.2f}  spec={r['specificity']:.2f}\n")

        # Worst performing
        f.write("\n── POOR PERFORMING (AUC < 0.55) " + "─" * 39 + "\n")
        for r in sorted(results, key=lambda x: x["auc"] if x["auc"] else 1):
            if r["auc"] and r["auc"] < 0.55:
                f.write(f"  {r['trait']:45s} AUC={r['auc']:.3f}  has={r['n_has']:3d}  no={r['n_no']:3d}  "
                        f"mean_diff={r['separation']:+5.1f}  {r['pgs_id']}  {r['flag']}\n")

        # Worse than random (candidates for blocklist)
        f.write("\n── BLOCKLIST CANDIDATES (AUC < 0.50) " + "─" * 34 + "\n")
        blocklist_candidates = [r for r in results if r["auc"] and r["auc"] < 0.50]
        if blocklist_candidates:
            for r in sorted(blocklist_candidates, key=lambda x: x["auc"]):
                f.write(f"  {r['trait']:45s} AUC={r['auc']:.3f}  {r['pgs_id']}  n_models={r['n_models']}\n")
        else:
            f.write("  None found\n")

        # Known problematic models
        f.write("\n── KNOWN PROBLEMATIC MODELS " + "─" * 43 + "\n")
        for trait_name in ["otosclerosis", "colorectal cancer"]:
            found = [r for r in results if r["trait"] == trait_name]
            if found:
                r = found[0]
                f.write(f"  {r['trait']:45s} AUC={r['auc']:.3f}  {r['pgs_id']}  n_models={r['n_models']}  "
                        f"mean_has={r['mean_has']}  mean_no={r['mean_no']}\n")
            else:
                f.write(f"  {trait_name:45s} (not enough data to validate)\n")

    print(f"  Summary written to: {summary_path}")

    # Print summary to stdout
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"Conditions validated: {len(results)}")
    if aucs:
        print(f"Mean AUC: {sum(aucs)/len(aucs):.3f}")
        print(f"Good (≥0.65): {sum(1 for a in aucs if a >= 0.65)}")
        print(f"Moderate (0.55-0.65): {sum(1 for a in aucs if 0.55 <= a < 0.65)}")
        print(f"Poor (<0.55): {sum(1 for a in aucs if a < 0.55)}")
        print(f"Worse than random (<0.50): {sum(1 for a in aucs if a < 0.50)}")

    print("\nTop 10 by AUC:")
    for r in results[:10]:
        print(f"  {r['trait']:40s} AUC={r['auc']:.3f}  n={r['n_has']}+{r['n_no']}  "
              f"sens={r['sensitivity']:.2f}  spec={r['specificity']:.2f}")

    print("\nWorst by AUC:")
    for r in sorted(results, key=lambda x: x["auc"] if x["auc"] else 1)[:5]:
        print(f"  {r['trait']:40s} AUC={r['auc']:.3f}  n={r['n_has']}+{r['n_no']}  {r['flag']}")

    blocklist = [r for r in results if r["auc"] and r["auc"] < 0.50]
    if blocklist:
        print(f"\nBLOCKLIST CANDIDATES ({len(blocklist)}):")
        for r in blocklist:
            print(f"  {r['trait']:40s} AUC={r['auc']:.3f}  {r['pgs_id']}")

    return results


if __name__ == "__main__":
    results = main()
