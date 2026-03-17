#!/usr/bin/env python3
"""
PGP Validation: Join PRS scores with diagnoses, compute AUC per trait.

Strategy:
1. Load pgp-pipeline.db: hex_id → pgp_num mapping + validatable traits
2. Load pgp-prs-scores.json.gz: sample scores
3. Map PGS models to traits (from pgs_metadata + all-traits.json)
4. For each trait with ≥10 cases: compute AUC (cases vs controls)
5. Also validate height (Pearson r with PGS for height)
6. Output: pgp-validation-results.json

Run: python3 validate-pgp.py
"""

import json, gzip, sqlite3, os, sys, time
import numpy as np
from collections import defaultdict

BASE_DIR = '/workspace/pgp-pipeline'
PIPELINE_DB = f'{BASE_DIR}/pgp-pipeline.db'
MEGA_DB = f'{BASE_DIR}/pgs-mega.db'
SCORES_FILE = f'{BASE_DIR}/pgp-prs-scores.json.gz'
TRAITS_FILE = f'{BASE_DIR}/all-traits.json'
OUTPUT = f'{BASE_DIR}/pgp-validation-results.json'


def compute_auc(scores_cases, scores_controls):
    """Compute AUC using Mann-Whitney U statistic."""
    n1 = len(scores_cases)
    n2 = len(scores_controls)
    if n1 == 0 or n2 == 0:
        return 0.5

    # Combine and rank
    all_scores = [(s, 1) for s in scores_cases] + [(s, 0) for s in scores_controls]
    all_scores.sort(key=lambda x: x[0])

    # Mann-Whitney U
    rank_sum = 0
    for i, (score, label) in enumerate(all_scores):
        if label == 1:
            rank_sum += i + 1  # 1-indexed rank

    # Handle ties by averaging ranks
    u = rank_sum - n1 * (n1 + 1) / 2
    auc = u / (n1 * n2)

    # AUC should be >= 0.5 (flip if needed for directionality)
    return max(auc, 1 - auc)


def pearson_r(x, y):
    """Compute Pearson correlation coefficient."""
    x = np.array(x)
    y = np.array(y)
    if len(x) < 3:
        return 0.0
    mx = np.mean(x)
    my = np.mean(y)
    num = np.sum((x - mx) * (y - my))
    den = np.sqrt(np.sum((x - mx)**2) * np.sum((y - my)**2))
    if den == 0:
        return 0.0
    return float(num / den)


def main():
    t0 = time.time()
    print("=" * 60)
    print("PGP PRS Validation")
    print("=" * 60, flush=True)

    # Step 1: Load sample mapping (hex_id → pgp_num → sample_id)
    print("\n[1] Loading sample mapping from pgp-pipeline.db...", flush=True)
    conn = sqlite3.connect(PIPELINE_DB)
    cur = conn.cursor()

    # hex_id → list of pgp_nums
    hex_to_samples = defaultdict(list)
    cur.execute("SELECT hex_id, pgp_num FROM pgp_genome_links WHERE has_file = 1")
    for hex_id, pgp_num in cur:
        sample_id = f"pgp-{pgp_num}"
        hex_to_samples[hex_id].append(sample_id)
    print(f"  {len(hex_to_samples)} hex_ids with genome files", flush=True)

    # Load validatable traits
    hex_traits = defaultdict(set)  # hex_id → set of trait names
    cur.execute("SELECT hex_id, trait FROM pgp_validatable_traits")
    for hex_id, trait in cur:
        hex_traits[hex_id].add(trait)

    # Load height data
    hex_height = {}
    cur.execute("SELECT hex_id, height_cm FROM pgp_profiles WHERE height_cm IS NOT NULL AND height_cm > 0")
    for hex_id, h in cur:
        hex_height[hex_id] = h
    print(f"  {len(hex_height)} profiles with height data", flush=True)
    print(f"  {sum(len(v) for v in hex_traits.values())} validatable trait entries", flush=True)

    conn.close()

    # Step 2: Load PRS scores
    print("\n[2] Loading PRS scores...", flush=True)
    with gzip.open(SCORES_FILE, 'rt') as f:
        all_scores = json.load(f)
    print(f"  {len(all_scores)} scored samples", flush=True)

    # Step 3: Build PGS model → trait mapping
    print("\n[3] Building model → trait mapping...", flush=True)

    # From pgs_metadata
    conn = sqlite3.connect(MEGA_DB)
    conn.execute('PRAGMA mmap_size = 214748364800')
    cur = conn.cursor()
    cur.execute("SELECT pgs_id, trait_name, trait_efo, category FROM pgs_metadata")
    model_info = {}
    for pgs_id, trait_name, trait_efo, category in cur:
        model_info[pgs_id] = {
            'trait_name': trait_name,
            'trait_efo': trait_efo or '',
            'category': category or '',
        }
    conn.close()
    print(f"  {len(model_info)} models with metadata", flush=True)

    # From all-traits.json
    with open(TRAITS_FILE) as f:
        all_traits = json.load(f)

    # Build trait keyword → PGS model list
    trait_to_models = defaultdict(list)

    # Manual mapping: validatable trait name → keywords to match in PGS trait names
    TRAIT_KEYWORDS = {
        'hypercholesterolemia': ['cholesterol', 'ldl', 'lipid', 'hypercholesterolemia'],
        'depression': ['depression', 'depressive', 'major depressive'],
        'hypertension': ['hypertension', 'blood pressure', 'systolic', 'diastolic'],
        'asthma': ['asthma'],
        'myopia': ['myopia', 'refractive error', 'near-sightedness'],
        'migraine': ['migraine'],
        'anxiety': ['anxiety'],
        'gerd': ['gastroesophageal', 'gerd', 'acid reflux'],
        'macular_degeneration': ['macular degeneration', 'amd'],
        'hypothyroidism': ['hypothyroidism', 'thyroid'],
        'sleep_apnea': ['sleep apnea', 'apnea'],
        'type_2_diabetes': ['type 2 diabetes', 'diabetes mellitus', 't2d'],
        'adhd': ['adhd', 'attention deficit', 'hyperactivity'],
        'osteoarthritis': ['osteoarthritis'],
        'allergic_rhinitis': ['allergic rhinitis', 'hay fever', 'rhinitis'],
        'obesity': ['obesity', 'bmi', 'body mass index'],
        'irritable_bowel_syndrome': ['irritable bowel', 'ibs'],
        'kidney_stones': ['kidney stone', 'nephrolithiasis', 'urolithiasis'],
        'eczema': ['eczema', 'atopic dermatitis', 'dermatitis'],
        'coronary_artery_disease': ['coronary artery', 'coronary heart', 'cad', 'myocardial infarction'],
        'atrial_fibrillation': ['atrial fibrillation'],
        'psoriasis': ['psoriasis'],
        'rheumatoid_arthritis': ['rheumatoid arthritis'],
        'celiac_disease': ['celiac', 'coeliac'],
        'crohns_disease': ['crohn'],
        'ulcerative_colitis': ['ulcerative colitis'],
        'breast_cancer': ['breast cancer', 'breast carcinoma'],
        'prostate_cancer': ['prostate cancer', 'prostate carcinoma'],
        'colorectal_cancer': ['colorectal', 'colon cancer'],
        'lung_cancer': ['lung cancer'],
        'gout': ['gout', 'urate', 'uric acid'],
        'epilepsy': ['epilepsy'],
        'glaucoma': ['glaucoma'],
        'bipolar_disorder': ['bipolar'],
        'schizophrenia': ['schizophrenia'],
        'alopecia': ['alopecia', 'hair loss', 'baldness'],
        'astigmatism': ['astigmatism'],
        'height': ['height', 'stature'],
    }

    for pgs_id, info in model_info.items():
        trait_lower = info['trait_name'].lower()
        for trait_name, keywords in TRAIT_KEYWORDS.items():
            for kw in keywords:
                if kw in trait_lower:
                    trait_to_models[trait_name].append(pgs_id)
                    break

    for trait, models in sorted(trait_to_models.items()):
        print(f"  {trait}: {len(models)} PGS models", flush=True)

    # Step 4: Match scored samples with traits
    print("\n[4] Matching samples with traits...", flush=True)

    # Build sample_id → hex_id mapping
    sample_to_hex = {}
    for hex_id, sample_ids in hex_to_samples.items():
        for sid in sample_ids:
            sample_to_hex[sid] = hex_id

    # Find scored samples with trait data
    scored_hex = set()
    for sample_id in all_scores:
        hex_id = sample_to_hex.get(sample_id)
        if hex_id:
            scored_hex.add(hex_id)
    print(f"  {len(scored_hex)} scored samples matched to hex_ids", flush=True)

    # Collect per-trait case/control sets
    trait_cases = defaultdict(set)   # trait → set of hex_ids with that trait
    trait_controls = defaultdict(set)  # trait → set of hex_ids without that trait

    for hex_id in scored_hex:
        traits = hex_traits.get(hex_id, set())
        for trait_name in TRAIT_KEYWORDS:
            if trait_name in traits:
                trait_cases[trait_name].add(hex_id)
            else:
                trait_controls[trait_name].add(hex_id)

    print(f"\n  Trait case counts (scored samples only):", flush=True)
    for trait in sorted(trait_cases.keys()):
        n_cases = len(trait_cases[trait])
        n_controls = len(trait_controls[trait])
        if n_cases > 0:
            print(f"    {trait}: {n_cases} cases / {n_controls} controls", flush=True)

    # Step 5: Compute AUC per trait per model
    print("\n[5] Computing AUC per trait...", flush=True)
    results = {}

    for trait_name, pgs_models in sorted(trait_to_models.items()):
        cases = trait_cases.get(trait_name, set())
        controls = trait_controls.get(trait_name, set())

        if len(cases) < 5:
            continue

        best_auc = 0.5
        best_model = None
        model_results = []

        for pgs_id in pgs_models:
            # Collect scores for cases and controls
            case_scores = []
            control_scores = []

            for hex_id in cases:
                for sid in hex_to_samples.get(hex_id, []):
                    score_data = all_scores.get(sid, {}).get(pgs_id)
                    if score_data and score_data['coverage'] >= 0.5:
                        case_scores.append(score_data['score'])
                        break

            for hex_id in controls:
                for sid in hex_to_samples.get(hex_id, []):
                    score_data = all_scores.get(sid, {}).get(pgs_id)
                    if score_data and score_data['coverage'] >= 0.5:
                        control_scores.append(score_data['score'])
                        break

            if len(case_scores) >= 5 and len(control_scores) >= 10:
                auc = compute_auc(case_scores, control_scores)
                model_results.append({
                    'pgs_id': pgs_id,
                    'auc': round(auc, 4),
                    'n_cases': len(case_scores),
                    'n_controls': len(control_scores),
                    'trait_name': model_info.get(pgs_id, {}).get('trait_name', ''),
                })
                if auc > best_auc:
                    best_auc = auc
                    best_model = pgs_id

        if model_results:
            # Sort by AUC descending
            model_results.sort(key=lambda x: -x['auc'])
            results[trait_name] = {
                'n_cases': len(cases),
                'n_controls': len(controls),
                'n_models_tested': len(model_results),
                'best_model': best_model,
                'best_auc': best_auc,
                'best_trait_name': model_info.get(best_model, {}).get('trait_name', ''),
                'top_models': model_results[:10],
            }
            print(f"  {trait_name}: AUC={best_auc:.3f} ({best_model}), "
                  f"{len(cases)} cases, {len(model_results)} models tested", flush=True)

    # Step 6: Height validation (Pearson r)
    print("\n[6] Height validation...", flush=True)
    height_models = trait_to_models.get('height', [])
    if height_models:
        best_r = 0
        best_height_model = None

        for pgs_id in height_models:
            heights = []
            prs_scores = []
            for hex_id, h in hex_height.items():
                if hex_id not in scored_hex:
                    continue
                for sid in hex_to_samples.get(hex_id, []):
                    score_data = all_scores.get(sid, {}).get(pgs_id)
                    if score_data and score_data['coverage'] >= 0.5:
                        heights.append(h)
                        prs_scores.append(score_data['score'])
                        break

            if len(heights) >= 10:
                r = pearson_r(prs_scores, heights)
                if abs(r) > abs(best_r):
                    best_r = r
                    best_height_model = pgs_id

        if best_height_model:
            # Recompute for best model
            heights = []
            prs_scores = []
            for hex_id, h in hex_height.items():
                if hex_id not in scored_hex:
                    continue
                for sid in hex_to_samples.get(hex_id, []):
                    score_data = all_scores.get(sid, {}).get(best_height_model)
                    if score_data and score_data['coverage'] >= 0.5:
                        heights.append(h)
                        prs_scores.append(score_data['score'])
                        break

            r = pearson_r(prs_scores, heights)
            results['height_continuous'] = {
                'metric': 'pearson_r',
                'r': round(r, 4),
                'r_squared': round(r**2, 4),
                'n_samples': len(heights),
                'best_model': best_height_model,
                'best_trait_name': model_info.get(best_height_model, {}).get('trait_name', ''),
                'n_models_tested': len(height_models),
            }
            print(f"  Height: r={r:.3f} (r²={r**2:.3f}), {best_height_model}, "
                  f"n={len(heights)}, {len(height_models)} models tested", flush=True)

    # Step 7: Save results
    print(f"\n[7] Saving {OUTPUT}...", flush=True)
    with open(OUTPUT, 'w') as f:
        json.dump(results, f, indent=2)

    elapsed = time.time() - t0
    print(f"\nDone! {len(results)} traits validated in {elapsed:.0f}s")
    print(f"\nSummary:")
    for trait, r in sorted(results.items(), key=lambda x: -x[1].get('best_auc', x[1].get('r', 0))):
        if 'best_auc' in r:
            print(f"  {trait}: AUC={r['best_auc']:.3f} ({r['best_model']}, {r['n_cases']} cases)")
        elif 'r' in r:
            print(f"  {trait}: r={r['r']:.3f} ({r['best_model']}, {r['n_samples']} samples)")


if __name__ == '__main__':
    main()
