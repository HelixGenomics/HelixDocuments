#!/usr/bin/env python3
"""
Combine PGP + OpenSNP validation results and select best PGS model per trait.
Output: best-models.json for production use.
"""

import json, os

BASE_DIR = '/workspace/pgp-pipeline'
PGP_RESULTS = f'{BASE_DIR}/pgp-validation-results.json'
OPENSNP_RESULTS = f'{BASE_DIR}/opensnp-validation-results-v2.json'
OUTPUT = f'{BASE_DIR}/best-models.json'


def main():
    print("=" * 60)
    print("Combining validation results — best model selection")
    print("=" * 60)

    # Load PGP results
    with open(PGP_RESULTS) as f:
        pgp = json.load(f)
    print(f"\nPGP: {len(pgp)} traits validated")

    # Load OpenSNP results
    with open(OPENSNP_RESULTS) as f:
        opensnp = json.load(f)

    # Collect all height models
    height_models = {}
    for h in opensnp.get('height', []):
        pid = h['pgs_id']
        height_models[pid] = {
            'r_opensnp': h.get('r_all', 0),
            'r_eur_opensnp': h.get('r_eur', 0),
            'n_opensnp': h.get('n_all', 0),
            'coverage_opensnp': h.get('coverage', 0),
        }

    # PGP height
    pgp_height = pgp.get('height_continuous', {})
    if pgp_height:
        pid = pgp_height['best_model']
        if pid not in height_models:
            height_models[pid] = {}
        height_models[pid]['r_pgp'] = pgp_height['r']
        height_models[pid]['n_pgp'] = pgp_height['n_samples']

    # Build best models
    best = {}

    # Height: combine OpenSNP and PGP results
    print("\n--- HEIGHT ---")
    best_height = None
    best_height_score = 0
    for pid, info in height_models.items():
        # Weighted combination: prefer OpenSNP (more samples) but consider PGP
        r_osn = abs(info.get('r_opensnp', 0))
        r_pgp = abs(info.get('r_pgp', 0))
        combined = r_osn * 0.7 + r_pgp * 0.3 if r_pgp else r_osn
        if combined > best_height_score:
            best_height_score = combined
            best_height = pid

    if best_height:
        info = height_models[best_height]
        best['height'] = {
            'pgs_id': best_height,
            'metric': 'pearson_r',
            'r_opensnp': info.get('r_opensnp', 0),
            'r_pgp': info.get('r_pgp', 0),
            'source': 'combined',
        }
        print(f"  Best: {best_height} (r_osn={info.get('r_opensnp',0):.3f}, r_pgp={info.get('r_pgp',0):.3f})")

    # Disease traits from PGP
    print("\n--- DISEASE TRAITS (PGP) ---")
    for trait, info in sorted(pgp.items()):
        if trait == 'height_continuous':
            continue
        if 'best_auc' not in info:
            continue
        auc = info['best_auc']
        pid = info['best_model']
        n_cases = info['n_cases']

        # Confidence: higher AUC with more cases = more confident
        confidence = 'high' if auc > 0.65 and n_cases >= 20 else \
                     'medium' if auc > 0.55 and n_cases >= 10 else 'low'

        best[trait] = {
            'pgs_id': pid,
            'metric': 'auc',
            'auc': auc,
            'n_cases': n_cases,
            'n_controls': info['n_controls'],
            'trait_name': info.get('best_trait_name', ''),
            'n_models_tested': info.get('n_models_tested', 0),
            'source': 'pgp',
            'confidence': confidence,
        }
        print(f"  {trait}: AUC={auc:.3f} ({pid}), {n_cases} cases [{confidence}]")

    # Hair color from OpenSNP
    print("\n--- HAIR COLOR (OpenSNP) ---")
    hair = opensnp.get('hair_color', [])
    if hair:
        # Group by target
        for target in ['black', 'brown', 'blonde', 'red']:
            target_models = [h for h in hair if h.get('target') == target]
            if target_models:
                best_hair = max(target_models, key=lambda x: x.get('auc', 0))
                trait_key = f'hair_{target}'
                best[trait_key] = {
                    'pgs_id': best_hair['pgs_id'],
                    'metric': 'auc',
                    'auc': best_hair['auc'],
                    'n_cases': best_hair.get('n_cases', 0),
                    'source': 'opensnp',
                    'trait_name': best_hair.get('trait', ''),
                }
                print(f"  hair_{target}: AUC={best_hair['auc']:.3f} ({best_hair['pgs_id']})")

    # Save
    print(f"\n\nTotal: {len(best)} traits with best models")
    with open(OUTPUT, 'w') as f:
        json.dump(best, f, indent=2)
    print(f"Saved to {OUTPUT}")

    # Summary table
    print("\n" + "=" * 70)
    print(f"{'Trait':<30} {'Model':<15} {'Metric':<10} {'Value':<8} {'Confidence':<10}")
    print("-" * 70)
    for trait in sorted(best.keys()):
        info = best[trait]
        pid = info['pgs_id']
        metric = info['metric']
        if metric == 'auc':
            val = f"{info['auc']:.3f}"
        elif metric == 'pearson_r':
            val = f"{max(abs(info.get('r_opensnp',0)), abs(info.get('r_pgp',0))):.3f}"
        else:
            val = "?"
        conf = info.get('confidence', 'n/a')
        print(f"  {trait:<28} {pid:<15} {metric:<10} {val:<8} {conf}")


if __name__ == '__main__':
    main()
