#!/usr/bin/env python3
"""Fetch ancestry metadata for all PGS models we have, build per-trait ancestry selection map."""

import os, json, time, gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime
from collections import defaultdict

COMPACT_DIR = '/workspace/pgp-pipeline/pgs-compact'
OUTPUT_META = '/workspace/1000g-scoring/pgs-metadata.json'
OUTPUT_SELECTION = '/workspace/1000g-scoring/pgs-ancestry-selection.json'

ANCESTRY_GROUPS = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR']

def ts():
    return datetime.now().strftime('%H:%M:%S')

def fetch_one(pgs_id):
    """Fetch metadata for a single PGS from the catalog API."""
    url = f'https://www.pgscatalog.org/rest/score/{pgs_id}'
    for attempt in range(3):
        try:
            req = Request(url, headers={'Accept': 'application/json'})
            with urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            
            # Extract what we need
            anc_dist = data.get('ancestry_distribution', {})
            gwas_dist = anc_dist.get('gwas', {}).get('dist', {})
            gwas_count = anc_dist.get('gwas', {}).get('count', 0)
            eval_dist = anc_dist.get('eval', {}).get('dist', {})
            dev_dist = anc_dist.get('dev', {}).get('dist', {})
            
            # Get trait EFO IDs
            trait_efo = []
            trait_reported = data.get('trait_reported', '')
            for t in data.get('trait_efo', []):
                trait_efo.append({
                    'id': t.get('id', ''),
                    'label': t.get('label', ''),
                })
            
            return {
                'pgs_id': pgs_id,
                'name': data.get('name', ''),
                'trait_reported': trait_reported,
                'trait_efo': trait_efo,
                'variants_number': data.get('variants_number', 0),
                'gwas_ancestry': gwas_dist,
                'gwas_sample_count': gwas_count,
                'eval_ancestry': eval_dist,
                'dev_ancestry': dev_dist,
                'publication': data.get('publication', {}).get('id', ''),
                'license': data.get('license', ''),
                # Compute diversity score: lower EUR% = more diverse
                'eur_gwas_pct': gwas_dist.get('EUR', 0),
                'multi_ancestry': len([k for k in gwas_dist if k in ('AFR','EAS','SAS','AMR','ASN') and gwas_dist[k] > 1]) > 0,
                'error': None,
            }
        except HTTPError as e:
            if e.code == 429:  # Rate limited
                time.sleep(5 * (attempt + 1))
                continue
            elif e.code == 404:
                return {'pgs_id': pgs_id, 'error': 'not_found'}
            else:
                time.sleep(2)
                continue
        except Exception as e:
            time.sleep(2)
            continue
    
    return {'pgs_id': pgs_id, 'error': 'failed_after_retries'}


def score_model(meta, target_ancestry):
    """Score a PGS model for a given target ancestry. Higher = better."""
    if meta.get('error'):
        return -1
    
    variants = meta.get('variants_number', 0)
    gwas = meta.get('gwas_ancestry', {})
    eval_anc = meta.get('eval_ancestry', {})
    gwas_n = meta.get('gwas_sample_count', 0)
    
    score = 0.0
    
    # 1. Variant count (log scale, max ~25 points)
    import math
    if variants > 0:
        score += min(25, math.log10(variants) * 4)
    
    # 2. Target ancestry in GWAS (max 30 points)
    # Map ancestry groups - some catalog entries use ASN instead of EAS
    target_pct = gwas.get(target_ancestry, 0)
    if target_ancestry == 'EAS':
        target_pct = max(target_pct, gwas.get('ASN', 0))
    if target_ancestry == 'ALL':
        # For ALL, reward diversity
        n_ancestries = len([k for k in gwas if gwas[k] > 1])
        target_pct = min(100, n_ancestries * 20)
    score += target_pct * 0.3
    
    # 3. Multi-ancestry bonus (max 15 points)
    n_diverse = len([k for k in gwas if k in ('AFR','EAS','SAS','AMR','ASN') and gwas[k] > 1])
    score += n_diverse * 3.75
    
    # 4. GWAS sample size (log scale, max 15 points)
    if gwas_n > 0:
        score += min(15, math.log10(gwas_n) * 2.5)
    
    # 5. Validated in target ancestry (max 15 points)
    if target_ancestry in eval_anc and eval_anc[target_ancestry] > 0:
        score += 15
    elif 'MAE' in eval_anc:  # Multi-ancestry evaluation
        score += 10
    
    return round(score, 2)


def build_selection_map(all_meta):
    """For each trait, select the best PGS per ancestry group."""
    
    # Group models by trait (EFO ID)
    trait_models = defaultdict(list)
    for m in all_meta:
        if m.get('error'):
            continue
        for t in m.get('trait_efo', []):
            efo_id = t['id']
            trait_models[efo_id].append(m)
        # Also group by reported trait name as fallback
        if m.get('trait_reported'):
            trait_models[f"reported:{m['trait_reported']}"].append(m)
    
    selection = {}
    
    for trait_key, models in trait_models.items():
        if len(models) == 0:
            continue
        
        trait_selection = {
            'trait_key': trait_key,
            'trait_label': models[0].get('trait_efo', [{}])[0].get('label', '') if 'reported:' not in trait_key else trait_key.replace('reported:', ''),
            'n_models': len(models),
            'best_per_ancestry': {},
            'all_models': [],
        }
        
        # Brief summary of each model for this trait
        for m in models:
            trait_selection['all_models'].append({
                'pgs_id': m['pgs_id'],
                'variants': m['variants_number'],
                'eur_pct': m['eur_gwas_pct'],
                'multi_ancestry': m['multi_ancestry'],
            })
        
        # Select best model per ancestry
        for anc in ANCESTRY_GROUPS + ['ALL']:
            scored = [(score_model(m, anc), m) for m in models]
            scored.sort(key=lambda x: -x[0])
            
            if scored and scored[0][0] > 0:
                best = scored[0][1]
                trait_selection['best_per_ancestry'][anc] = {
                    'pgs_id': best['pgs_id'],
                    'score': scored[0][0],
                    'variants': best['variants_number'],
                    'gwas_ancestry': best['gwas_ancestry'],
                }
                # If runner-up is significantly different, note it
                if len(scored) > 1 and scored[1][0] > scored[0][0] * 0.8:
                    trait_selection['best_per_ancestry'][anc]['runner_up'] = {
                        'pgs_id': scored[1][1]['pgs_id'],
                        'score': scored[1][0],
                    }
        
        selection[trait_key] = trait_selection
    
    return selection


def main():
    print(f'[{ts()}] PGS Ancestry Metadata Fetcher', flush=True)
    
    # Get list of PGS IDs from compact files
    pgs_ids = sorted([f.replace('.gz', '') for f in os.listdir(COMPACT_DIR) if f.endswith('.gz')])
    print(f'  {len(pgs_ids)} PGS models to fetch metadata for', flush=True)
    
    # Parallel fetch with rate limiting (10 workers to be nice to API)
    all_meta = []
    done = 0
    errors = 0
    
    with ThreadPoolExecutor(max_workers=10) as pool:
        futures = {pool.submit(fetch_one, pid): pid for pid in pgs_ids}
        for future in as_completed(futures):
            result = future.result()
            all_meta.append(result)
            done += 1
            if result.get('error'):
                errors += 1
            
            if done % 100 == 0:
                print(f'  [{ts()}] Fetched {done}/{len(pgs_ids)} ({errors} errors)', flush=True)
            
            # Small delay to avoid rate limiting
            if done % 50 == 0:
                time.sleep(1)
    
    # Sort by PGS ID
    all_meta.sort(key=lambda m: m.get('pgs_id', ''))
    
    # Save full metadata
    os.makedirs(os.path.dirname(OUTPUT_META), exist_ok=True)
    with open(OUTPUT_META, 'w') as f:
        json.dump(all_meta, f, indent=2)
    print(f'\n  Full metadata saved to {OUTPUT_META}', flush=True)
    
    # Summary stats
    valid = [m for m in all_meta if not m.get('error')]
    multi = [m for m in valid if m.get('multi_ancestry')]
    eur_only = [m for m in valid if m.get('eur_gwas_pct', 0) >= 95]
    no_gwas = [m for m in valid if not m.get('gwas_ancestry')]
    
    print(f'\n{"="*60}', flush=True)
    print(f'PGS Ancestry Metadata Summary — {len(valid)} models', flush=True)
    print(f'{"="*60}', flush=True)
    print(f'  Multi-ancestry GWAS (>1% non-EUR): {len(multi)}', flush=True)
    print(f'  EUR-dominant (>=95% EUR GWAS): {len(eur_only)}', flush=True)
    print(f'  No GWAS ancestry data: {len(no_gwas)}', flush=True)
    print(f'  Fetch errors: {errors}', flush=True)
    
    # Variant count distribution
    var_counts = [m['variants_number'] for m in valid if m.get('variants_number')]
    if var_counts:
        var_counts.sort()
        print(f'\n  Variant counts:', flush=True)
        print(f'    Min: {min(var_counts):,}', flush=True)
        print(f'    Median: {var_counts[len(var_counts)//2]:,}', flush=True)
        print(f'    Mean: {sum(var_counts)//len(var_counts):,}', flush=True)
        print(f'    Max: {max(var_counts):,}', flush=True)
        print(f'    >100K variants: {sum(1 for v in var_counts if v > 100000)}', flush=True)
        print(f'    >1M variants: {sum(1 for v in var_counts if v > 1000000)}', flush=True)
    
    # Build ancestry selection map
    print(f'\n[{ts()}] Building ancestry selection map...', flush=True)
    selection = build_selection_map(valid)
    
    # Count traits with different best models per ancestry
    traits_with_diff = 0
    for trait_key, sel in selection.items():
        if trait_key.startswith('reported:'):
            continue  # Skip duplicate reported-name entries
        best_ids = set(v['pgs_id'] for v in sel['best_per_ancestry'].values())
        if len(best_ids) > 1:
            traits_with_diff += 1
    
    efo_traits = {k: v for k, v in selection.items() if not k.startswith('reported:')}
    print(f'  {len(efo_traits)} unique EFO traits', flush=True)
    print(f'  {traits_with_diff} traits where best model differs by ancestry', flush=True)
    
    # Show examples where ancestry matters
    print(f'\n  Examples where ancestry changes best model:', flush=True)
    count = 0
    for trait_key, sel in sorted(selection.items()):
        if trait_key.startswith('reported:'):
            continue
        bpa = sel['best_per_ancestry']
        if len(set(v['pgs_id'] for v in bpa.values())) > 1:
            eur_best = bpa.get('EUR', {}).get('pgs_id', '?')
            afr_best = bpa.get('AFR', {}).get('pgs_id', '?')
            eas_best = bpa.get('EAS', {}).get('pgs_id', '?')
            label = sel.get('trait_label', trait_key)
            if eur_best != afr_best or eur_best != eas_best:
                print(f'    {label}: EUR→{eur_best}, AFR→{afr_best}, EAS→{eas_best}', flush=True)
                count += 1
                if count >= 15:
                    break
    
    with open(OUTPUT_SELECTION, 'w') as f:
        json.dump(selection, f, indent=2)
    print(f'\n  Selection map saved to {OUTPUT_SELECTION}', flush=True)
    
    print(f'\n[{ts()}] Done.', flush=True)


if __name__ == '__main__':
    main()
