#!/usr/bin/env python3
"""Scan full PGS Catalog for multi-ancestry models we don't have, then download them."""

import os, json, time, gzip
from urllib.request import urlopen, Request
from urllib.error import HTTPError
from datetime import datetime
from collections import defaultdict

COMPACT_DIR = '/workspace/pgp-pipeline/pgs-compact'
OUTPUT = '/workspace/1000g-scoring/missing-ancestry-models.json'
DOWNLOAD_LIST = '/workspace/1000g-scoring/ancestry-download-list.json'

def ts():
    return datetime.now().strftime('%H:%M:%S')

def main():
    # What we already have
    our_ids = set(f.replace('.gz','') for f in os.listdir(COMPACT_DIR) if f.endswith('.gz'))
    print(f'[{ts()}] We have {len(our_ids)} compact files', flush=True)
    
    # Scan full catalog
    all_scores = []
    offset = 0
    limit = 100
    total = None
    
    while True:
        url = f'https://www.pgscatalog.org/rest/score/all?limit={limit}&offset={offset}'
        req = Request(url, headers={'Accept': 'application/json'})
        
        for attempt in range(3):
            try:
                resp = urlopen(req, timeout=30)
                data = json.loads(resp.read().decode())
                break
            except HTTPError as e:
                if e.code == 429:
                    time.sleep(10)
                    continue
                raise
            except:
                time.sleep(3)
                continue
        else:
            print(f'  Failed at offset {offset}, stopping', flush=True)
            break
        
        if total is None:
            total = data.get('count', 0)
            print(f'  Total in catalog: {total}', flush=True)
        
        results = data.get('results', [])
        if not results:
            break
        
        for s in results:
            pgs_id = s['id']
            anc = s.get('ancestry_distribution', {})
            gwas = anc.get('gwas', {}).get('dist', {})
            dev = anc.get('dev', {}).get('dist', {})
            eval_d = anc.get('eval', {}).get('dist', {})
            gwas_n = anc.get('gwas', {}).get('count', 0)
            
            # Combine all ancestry sources
            all_anc = {}
            for d in [gwas, dev, eval_d]:
                for k, v in d.items():
                    if k not in all_anc or v > all_anc[k]:
                        all_anc[k] = v
            
            # Check for non-EUR ancestry presence
            non_eur_keys = ['AFR', 'EAS', 'SAS', 'AMR', 'ASN', 'MAE', 'MAO']
            has_non_eur = any(all_anc.get(k, 0) > 5 for k in non_eur_keys)
            
            # Check if developed/evaluated in non-EUR
            dev_non_eur = any(dev.get(k, 0) > 5 for k in non_eur_keys)
            eval_non_eur = any(eval_d.get(k, 0) > 5 for k in non_eur_keys)
            gwas_non_eur = any(gwas.get(k, 0) > 5 for k in non_eur_keys)
            
            trait_efo = [t.get('id','') for t in s.get('trait_efo', [])]
            
            all_scores.append({
                'id': pgs_id,
                'name': s.get('name', ''),
                'trait': s.get('trait_reported', ''),
                'trait_efo': trait_efo,
                'variants': s.get('variants_number', 0),
                'gwas_ancestry': gwas,
                'dev_ancestry': dev,
                'eval_ancestry': eval_d,
                'gwas_n': gwas_n,
                'has_non_eur': has_non_eur,
                'gwas_non_eur': gwas_non_eur,
                'dev_non_eur': dev_non_eur,
                'eval_non_eur': eval_non_eur,
                'we_have': pgs_id in our_ids,
                'ftp_url': s.get('ftp_scoring_file', ''),
            })
        
        offset += len(results)
        if offset % 500 == 0:
            print(f'  [{ts()}] Scanned {offset}/{total}...', flush=True)
        time.sleep(0.3)
    
    print(f'\n[{ts()}] Scanned {len(all_scores)} total scores', flush=True)
    
    # Analysis
    we_have = [s for s in all_scores if s['we_have']]
    missing = [s for s in all_scores if not s['we_have']]
    
    missing_with_ancestry = [s for s in missing if s['has_non_eur']]
    missing_gwas_diverse = [s for s in missing if s['gwas_non_eur']]
    missing_dev_diverse = [s for s in missing if s['dev_non_eur']]
    missing_eval_diverse = [s for s in missing if s['eval_non_eur']]
    
    print(f'\n{"="*60}', flush=True)
    print(f'CATALOG ANALYSIS', flush=True)
    print(f'{"="*60}', flush=True)
    print(f'Total in catalog: {len(all_scores)}', flush=True)
    print(f'We have: {len(we_have)}', flush=True)
    print(f'Missing: {len(missing)}', flush=True)
    print(f'', flush=True)
    print(f'Missing with ANY non-EUR ancestry: {len(missing_with_ancestry)}', flush=True)
    print(f'  - GWAS non-EUR: {len(missing_gwas_diverse)}', flush=True)
    print(f'  - Dev non-EUR: {len(missing_dev_diverse)}', flush=True)
    print(f'  - Eval non-EUR: {len(missing_eval_diverse)}', flush=True)
    
    # What ancestries are represented in missing models?
    anc_counts = defaultdict(int)
    for s in missing_with_ancestry:
        for src in [s['gwas_ancestry'], s['dev_ancestry'], s['eval_ancestry']]:
            for k in src:
                if src[k] > 5:
                    anc_counts[k] += 1
    
    print(f'\nAncestry representation in missing models:', flush=True)
    for k, v in sorted(anc_counts.items(), key=lambda x: -x[1]):
        print(f'  {k}: {v} models', flush=True)
    
    # What traits do the missing ancestry models cover?
    trait_counts = defaultdict(int)
    for s in missing_with_ancestry:
        trait_counts[s['trait'][:50]] += 1
    
    print(f'\nTop traits with missing ancestry models:', flush=True)
    for trait, count in sorted(trait_counts.items(), key=lambda x: -x[1])[:25]:
        print(f'  {trait}: {count}', flush=True)
    
    # Build download list — all missing models that have non-EUR ancestry
    to_download = []
    for s in missing_with_ancestry:
        if s['ftp_url']:
            to_download.append({
                'id': s['id'],
                'trait': s['trait'][:60],
                'variants': s['variants'],
                'url': s['ftp_url'],
                'gwas_ancestry': s['gwas_ancestry'],
                'dev_ancestry': s['dev_ancestry'],
            })
    
    to_download.sort(key=lambda x: x['id'])
    
    print(f'\nModels to download: {len(to_download)}', flush=True)
    
    with open(OUTPUT, 'w') as f:
        json.dump({'all_scores': all_scores}, f)
    
    with open(DOWNLOAD_LIST, 'w') as f:
        json.dump(to_download, f, indent=2)
    
    print(f'Download list saved to {DOWNLOAD_LIST}', flush=True)
    print(f'[{ts()}] Done.', flush=True)

if __name__ == '__main__':
    main()
