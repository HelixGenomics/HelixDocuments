#!/usr/bin/env python3
"""Download and convert missing ancestry PGS models to compact format."""

import os, json, gzip, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import urlopen, Request
from urllib.error import HTTPError
from datetime import datetime
from io import BytesIO, TextIOWrapper

DOWNLOAD_LIST = '/workspace/1000g-scoring/ancestry-download-list.json'
COMPACT_DIR = '/workspace/pgp-pipeline/pgs-compact'

def ts():
    return datetime.now().strftime('%H:%M:%S')

def download_and_convert(entry):
    """Download a PGS scoring file and convert to compact format."""
    pgs_id = entry['id']
    url = entry['url']
    out_path = os.path.join(COMPACT_DIR, f'{pgs_id}.gz')
    
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return pgs_id, 'skip', 0
    
    try:
        req = Request(url, headers={'User-Agent': 'HelixGenomics/1.0'})
        resp = urlopen(req, timeout=60)
        raw = resp.read()
        
        # Parse the scoring file (gzipped TSV from FTP)
        rows = []
        header_map = {}
        
        with gzip.open(BytesIO(raw), 'rt', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                
                # Header line
                if not header_map:
                    for i, col in enumerate(parts):
                        header_map[col.lower()] = i
                    continue
                
                # Need: rsid (or chr_name+chr_position), effect_allele, other_allele, effect_weight
                rsid = None
                ea = None
                oa = None
                weight = None
                
                # Get rsID
                if 'rsid' in header_map:
                    idx = header_map['rsid']
                    if idx < len(parts):
                        val = parts[idx].strip()
                        if val.startswith('rs'):
                            rsid = val
                
                # Get chr:pos as fallback
                if not rsid and 'chr_name' in header_map and 'chr_position' in header_map:
                    chr_idx = header_map['chr_name']
                    pos_idx = header_map['chr_position']
                    if chr_idx < len(parts) and pos_idx < len(parts):
                        chrom = parts[chr_idx].strip()
                        pos = parts[pos_idx].strip()
                        if chrom and pos:
                            rsid = f'chr{chrom}:{pos}'
                
                if not rsid:
                    continue
                
                # Effect allele
                if 'effect_allele' in header_map:
                    idx = header_map['effect_allele']
                    if idx < len(parts):
                        ea = parts[idx].strip().upper()
                
                # Other allele
                if 'other_allele' in header_map:
                    idx = header_map['other_allele']
                    if idx < len(parts):
                        oa = parts[idx].strip().upper()
                
                # Weight
                if 'effect_weight' in header_map:
                    idx = header_map['effect_weight']
                    if idx < len(parts):
                        try:
                            weight = float(parts[idx].strip())
                        except:
                            continue
                
                if not ea or weight is None:
                    continue
                
                if ea not in ('A', 'T', 'C', 'G'):
                    continue
                
                oa_clean = oa if oa and oa in ('A', 'T', 'C', 'G') else ''
                rows.append(f'{rsid}\t{ea}\t{oa_clean}\t{weight}')
        
        if not rows:
            return pgs_id, 'empty', 0
        
        # Write compact format
        with gzip.open(out_path, 'wt') as f:
            f.write('\n'.join(rows) + '\n')
        
        return pgs_id, 'ok', len(rows)
    
    except HTTPError as e:
        return pgs_id, f'http_{e.code}', 0
    except Exception as e:
        return pgs_id, f'error:{str(e)[:50]}', 0


def main():
    print(f'[{ts()}] Downloading missing ancestry PGS models', flush=True)
    
    with open(DOWNLOAD_LIST) as f:
        to_download = json.load(f)
    
    # Skip ones we already have
    existing = set(f.replace('.gz','') for f in os.listdir(COMPACT_DIR) if f.endswith('.gz'))
    to_download = [d for d in to_download if d['id'] not in existing]
    
    print(f'  {len(to_download)} models to download (after skipping existing)', flush=True)
    
    done = 0
    ok = 0
    fail = 0
    empty = 0
    total_variants = 0
    
    # 20 parallel workers
    with ThreadPoolExecutor(max_workers=20) as pool:
        futures = {pool.submit(download_and_convert, d): d for d in to_download}
        for future in as_completed(futures):
            pgs_id, status, n_rows = future.result()
            done += 1
            
            if status == 'ok':
                ok += 1
                total_variants += n_rows
            elif status == 'skip':
                ok += 1
            elif status == 'empty':
                empty += 1
            else:
                fail += 1
            
            if done % 100 == 0:
                print(f'  [{ts()}] {done}/{len(to_download)} ({ok} ok, {fail} fail, {empty} empty)', flush=True)
    
    # Final count
    final_count = len([f for f in os.listdir(COMPACT_DIR) if f.endswith('.gz')])
    print(f'\n{"="*60}', flush=True)
    print(f'Download complete', flush=True)
    print(f'  Downloaded: {ok} ok, {fail} failed, {empty} empty', flush=True)
    print(f'  Total variants added: {total_variants:,}', flush=True)
    print(f'  Total compact files now: {final_count}', flush=True)
    print(f'[{ts()}] Done.', flush=True)

if __name__ == '__main__':
    main()
