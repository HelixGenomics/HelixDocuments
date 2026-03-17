#!/usr/bin/env python3
"""
Phase 2: Score 1000G samples against PGS models from mega DB.
Produces per-ancestry population distributions for all needed models.

Strategy:
1. Load pgs-ancestry-selection.json to get the ~1,900 needed PGS IDs
2. For each chromosome, stream the 1000G VCF and build a position→dosage lookup
3. For each PGS model, query its variants from mega DB by chr:pos, look up dosages
4. Output: prs-distributions-mega.json.gz with per-ancestry sorted scores

Run on build machine: python3 score-1kg-mega.py [--workers N] [--chr N]
"""

import sqlite3
import json
import gzip
import os
import sys
import time
import struct
import mmap
import argparse
from collections import defaultdict
from multiprocessing import Pool, Manager, cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess

MEGA_DB = '/workspace/pgp-pipeline/pgs-mega.db'
ANCESTRY_SEL = '/workspace/1000g-scoring/pgs-ancestry-selection.json'
VCF_DIR = '/workspace/1000g-scoring/1000g'
PANEL_FILE = '/workspace/1000g-scoring/1000g/panel.txt'
OUTPUT = '/workspace/1000g-scoring/prs-distributions-mega.json.gz'
DOSES_DIR = '/workspace/1000g-scoring/1kg-doses'
POPS = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def load_sample_info():
    """Load 1000G sample names from VCF header + population from panel."""
    # Get sample names from a VCF header
    vcf_path = os.path.join(VCF_DIR, 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')
    result = subprocess.run(
        f"zcat {vcf_path} | head -300 | grep '^#CHROM'",
        shell=True, capture_output=True, text=True
    )
    header = result.stdout.strip().split('\t')
    sample_names = header[9:]
    n_samples = len(sample_names)
    print(f"  {n_samples} samples from VCF header")

    # Load population panel
    sample_pop = {}
    with open(PANEL_FILE) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[0] != 'sample':
                sample_pop[parts[0]] = parts[2]  # super_pop

    # Build per-population indices
    pop_indices = {pop: [] for pop in POPS}
    for i, name in enumerate(sample_names):
        pop = sample_pop.get(name)
        if pop and pop in pop_indices:
            pop_indices[pop].append(i)

    for pop in POPS:
        print(f"    {pop}: {len(pop_indices[pop])} samples")

    return sample_names, pop_indices, n_samples


def get_needed_pgs_ids():
    """Get the set of PGS IDs we need from ancestry selection."""
    data = json.load(open(ANCESTRY_SEL))
    needed = set()
    for trait_key, info in data.items():
        if not isinstance(info, dict) or 'best_per_ancestry' not in info:
            continue
        for pop, model in info['best_per_ancestry'].items():
            needed.add(model['pgs_id'])
            if 'runner_up' in model and isinstance(model['runner_up'], dict):
                rid = model['runner_up'].get('pgs_id')
                if rid:
                    needed.add(rid)
    return needed


def load_model_variants(pgs_id, db_path=MEGA_DB):
    """Load all variants for a PGS model from mega DB."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("""
        SELECT chr_name, chr_position, effect_allele, other_allele, effect_weight
        FROM pgs_variants WHERE pgs_id = ?
    """, (pgs_id,))
    variants = cur.fetchall()
    conn.close()
    return variants


def extract_chr_doses(chr_num, n_samples, needed_positions=None):
    """
    Stream a 1000G VCF for one chromosome, extract GT dosages.
    Returns: dict of (pos) -> [dosages for all samples]
    Only extracts positions in needed_positions if provided.
    """
    vcf_path = os.path.join(VCF_DIR, f'ALL.chr{chr_num}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')
    if not os.path.exists(vcf_path):
        print(f"  WARNING: VCF not found for chr{chr_num}")
        return {}

    doses = {}  # pos -> (ref, alt, dosages)
    count = 0
    skipped = 0

    proc = subprocess.Popen(
        ['zcat', vcf_path],
        stdout=subprocess.PIPE,
        bufsize=1024*1024*16
    )

    for raw_line in proc.stdout:
        line = raw_line.decode('utf-8', errors='replace')
        if line.startswith('#'):
            continue

        parts = line.split('\t', 10)  # only split first 10 fields
        if len(parts) < 10:
            continue

        pos = parts[1]

        # Filter to needed positions if set provided
        if needed_positions is not None and pos not in needed_positions:
            skipped += 1
            continue

        ref = parts[3].upper()
        alt = parts[4].upper()

        # Skip multiallelic
        if ',' in alt:
            skipped += 1
            continue

        # Parse GT field (format: 0|0, 0|1, 1|0, 1|1)
        gt_data = parts[9].rstrip('\n')
        sample_fields = gt_data.split('\t')

        dosages = []
        for sf in sample_fields:
            gt = sf.split(':')[0] if ':' in sf else sf
            if '|' in gt:
                a, b = gt.split('|')
            elif '/' in gt:
                a, b = gt.split('/')
            else:
                dosages.append(0.0)
                continue
            try:
                dosages.append(float(a) + float(b))
            except ValueError:
                dosages.append(0.0)

        if len(dosages) == n_samples:
            doses[pos] = (ref, alt, dosages)
            count += 1

        if count % 100000 == 0 and count > 0:
            print(f"    chr{chr_num}: {count} variants extracted, {skipped} skipped")

    proc.wait()
    print(f"  chr{chr_num}: {count} variants extracted total")
    return doses


def get_needed_positions_per_chr(needed_pgs_ids):
    """Query mega DB to find all needed chr:pos pairs per chromosome."""
    conn = sqlite3.connect(MEGA_DB)
    cur = conn.cursor()

    # Get unique positions needed per chromosome
    chr_positions = defaultdict(set)  # chr -> set of positions
    pgs_variants = defaultdict(list)  # pgs_id -> [(chr, pos, ea, oa, weight)]

    needed_list = sorted(needed_pgs_ids)
    total = len(needed_list)

    for i, pgs_id in enumerate(needed_list):
        cur.execute("""
            SELECT chr_name, chr_position, effect_allele, other_allele, effect_weight
            FROM pgs_variants WHERE pgs_id = ?
        """, (pgs_id,))

        for row in cur.fetchall():
            chr_name, pos, ea, oa, weight = row
            if chr_name and pos:
                # Normalize chr name (remove 'chr' prefix if present)
                c = chr_name.replace('chr', '')
                chr_positions[c].add(pos)
                pgs_variants[pgs_id].append((c, pos, ea or '', oa or '', weight))

        if (i + 1) % 200 == 0:
            print(f"  Loaded variants for {i+1}/{total} models...")

    conn.close()

    total_positions = sum(len(v) for v in chr_positions.values())
    print(f"\nTotal unique positions needed: {total_positions}")
    for c in sorted(chr_positions.keys(), key=lambda x: int(x) if x.isdigit() else 99):
        print(f"  chr{c}: {len(chr_positions[c])} positions")

    return chr_positions, pgs_variants


def score_model(pgs_id, model_variants, chr_doses, n_samples):
    """
    Score all samples for one PGS model.
    model_variants: [(chr, pos, ea, oa, weight)]
    chr_doses: {chr: {pos: (ref, alt, dosages)}}
    Returns: list of scores (length n_samples), matched count
    """
    scores = [0.0] * n_samples
    matched = 0
    total = len(model_variants)

    for chr_name, pos, ea, oa, weight in model_variants:
        dose_data = chr_doses.get(chr_name, {}).get(pos)
        if dose_data is None:
            continue

        ref, alt, dosages = dose_data

        # Allele alignment
        ea_upper = ea.upper()
        oa_upper = oa.upper()
        ref_upper = ref.upper()
        alt_upper = alt.upper()

        flip = False
        if ea_upper == alt_upper:
            flip = False  # dosage is already for ALT allele
        elif ea_upper == ref_upper:
            flip = True   # need to flip: use (2 - dosage)
        elif COMPLEMENT.get(ea_upper) == alt_upper:
            flip = False  # strand flip, ALT matches complement
        elif COMPLEMENT.get(ea_upper) == ref_upper:
            flip = True   # strand flip, REF matches complement
        else:
            continue  # can't align

        matched += 1
        for j in range(n_samples):
            d = dosages[j]
            if flip:
                d = 2.0 - d
            scores[j] += d * weight

    return scores, matched, total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=min(22, cpu_count()))
    parser.add_argument('--chr', type=str, help='Score single chromosome (for testing)')
    parser.add_argument('--resume', action='store_true', help='Resume from saved doses')
    args = parser.parse_args()

    t0 = time.time()
    print("=" * 60)
    print("1000G × Mega DB PRS Scoring Pipeline")
    print("=" * 60)

    # Step 1: Load sample info
    print("\n[Step 1] Loading sample info...")
    sample_names, pop_indices, n_samples = load_sample_info()

    # Step 2: Get needed PGS IDs
    print("\n[Step 2] Getting needed PGS IDs...")
    needed_pgs = get_needed_pgs_ids()
    # Also check which are actually in mega DB
    conn = sqlite3.connect(MEGA_DB)
    cur = conn.cursor()
    available = set()
    for pgs_id in needed_pgs:
        cur.execute("SELECT 1 FROM pgs_meta WHERE pgs_id = ?", (pgs_id,))
        if cur.fetchone():
            available.add(pgs_id)
    conn.close()
    needed_pgs = available
    print(f"  {len(needed_pgs)} PGS models needed and available in mega DB")

    # Step 3: Get needed positions per chromosome
    print("\n[Step 3] Loading model variants and needed positions...")
    chr_positions, pgs_variants = get_needed_positions_per_chr(needed_pgs)

    # Filter to valid chromosomes (1-22)
    valid_chrs = [str(i) for i in range(1, 23)]
    if args.chr:
        valid_chrs = [args.chr]

    # Step 4: Extract 1000G dosages per chromosome
    print(f"\n[Step 4] Extracting 1000G dosages ({len(valid_chrs)} chromosomes)...")
    os.makedirs(DOSES_DIR, exist_ok=True)

    chr_doses = {}  # chr -> {pos: (ref, alt, dosages)}

    for c in valid_chrs:
        cache_file = os.path.join(DOSES_DIR, f'chr{c}-doses.json.gz')

        if args.resume and os.path.exists(cache_file):
            print(f"  Loading cached chr{c} doses...")
            with gzip.open(cache_file, 'rt') as f:
                chr_doses[c] = json.load(f)
            # Convert dosage lists back
            for pos in chr_doses[c]:
                d = chr_doses[c][pos]
                chr_doses[c][pos] = (d[0], d[1], d[2])
            print(f"    chr{c}: {len(chr_doses[c])} variants from cache")
            continue

        needed = chr_positions.get(c, set())
        if not needed:
            print(f"  chr{c}: no positions needed, skipping")
            continue

        print(f"  Extracting chr{c} ({len(needed)} positions needed)...")
        doses = extract_chr_doses(c, n_samples, needed)
        chr_doses[c] = doses

        # Cache to disk
        print(f"    Caching chr{c} doses ({len(doses)} variants)...")
        with gzip.open(cache_file, 'wt') as f:
            # Convert to serializable format
            serializable = {}
            for pos, (ref, alt, dosages) in doses.items():
                serializable[pos] = [ref, alt, dosages]
            json.dump(serializable, f)

    total_variants = sum(len(d) for d in chr_doses.values())
    print(f"\n  Total 1000G variants loaded: {total_variants}")

    # Step 5: Score all models
    print(f"\n[Step 5] Scoring {len(needed_pgs)} models × {n_samples} samples...")
    distributions = {}
    scored = 0
    low_coverage = 0

    for pgs_id in sorted(needed_pgs):
        variants = pgs_variants.get(pgs_id, [])
        if not variants:
            continue

        scores, matched, total = score_model(pgs_id, variants, chr_doses, n_samples)

        if matched == 0:
            continue

        coverage = matched / total if total > 0 else 0

        if coverage < 0.05:
            low_coverage += 1
            continue

        # Build per-ancestry distributions
        dist = {}
        all_scores = sorted(scores)
        dist['ALL'] = {
            'scores': [round(s, 8) for s in all_scores],
            'n': n_samples,
            'mean': round(sum(scores) / n_samples, 8),
            'std': round((sum((s - sum(scores)/n_samples)**2 for s in scores) / n_samples) ** 0.5, 8),
            'matched': matched,
            'total': total,
            'coverage': round(coverage, 4),
        }

        for pop in POPS:
            idx = pop_indices[pop]
            if not idx:
                continue
            pop_scores = sorted([scores[i] for i in idx])
            n = len(pop_scores)
            mean = sum(pop_scores) / n
            std = (sum((s - mean)**2 for s in pop_scores) / n) ** 0.5
            dist[pop] = {
                'scores': [round(s, 8) for s in pop_scores],
                'n': n,
                'mean': round(mean, 8),
                'std': round(std, 8),
            }

        distributions[pgs_id] = dist
        scored += 1

        if scored % 100 == 0:
            elapsed = time.time() - t0
            print(f"  Scored {scored}/{len(needed_pgs)} models ({elapsed:.0f}s elapsed)")

    print(f"\n  Scored: {scored} models")
    print(f"  Low coverage (skipped): {low_coverage}")

    # Step 6: Write output
    print(f"\n[Step 6] Writing {OUTPUT}...")
    with gzip.open(OUTPUT, 'wt') as f:
        json.dump(distributions, f)

    size_mb = os.path.getsize(OUTPUT) / 1024 / 1024
    elapsed = time.time() - t0
    print(f"\nDone! {size_mb:.1f} MB written in {elapsed:.0f}s")
    print(f"Models with distributions: {len(distributions)}")


if __name__ == '__main__':
    main()
