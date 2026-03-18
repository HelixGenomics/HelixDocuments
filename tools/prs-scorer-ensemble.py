#!/usr/bin/env python3
"""
PRS scorer v6 — accurate scoring with EUR freq-based ambiguous SNP resolution + ensemble.

Key improvements over v5:
- Float64 dosages (was float32 — accumulated precision errors)
- Ambiguous A/T and C/G SNPs resolved using EUR allele frequencies
  (without this, ~8% of SNPs randomly scored in wrong direction)
- Ensemble scoring: groups models by trait, removes outliers, weighted z-average
- 32 parallel workers
"""
import os, sys, json, gzip, time, bisect, statistics, math, csv, re
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import sqlite3
import multiprocessing

SCORE_DIR = "/opt/helix/data/pgs-score-files"
TRAIT_FILE = "/opt/helix/data/plink2-trait-info.json"
DIST_FILE = "/opt/helix/data/prs-distributions-matched.json.gz"
BATCH_SIZE = 400
METADATA_CSV = "/opt/helix/data/pgs_all_metadata_scores.csv"
FREQ_STATS_CACHE = "/opt/helix/data/pgs-freq-stats-cache.json"
EUR_FREQ_FILE = "/opt/helix/data/prs-eur-freqs-1kg.json.gz"

BASE_ENCODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
COMP_TABLE = np.array([3, 2, 1, 0, 255], dtype=np.uint8)

# Shared arrays (set before fork, read-only in children via COW)
_RSID_ARR = None    # int64, sorted — genotype rsids
_DOS_ARR = None     # float64 — genotype dosages (was float32 in v5!)
_REF_ARR = None     # uint8 — genotype ref allele codes
_ALT_ARR = None     # uint8 — genotype alt allele codes
_EUR_RSID = None    # int64, sorted — EUR freq rsids
_EUR_ALLELE = None  # uint8 — EUR freq allele codes
_EUR_FREQ = None    # float32 — EUR frequencies


def load_genotypes_numpy(db_path):
    """Load genotypes from SQLite into sorted numpy arrays. Uses float64 for dosages."""
    t0 = time.time()
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA mmap_size = 1073741824")

    n_total = db.execute(
        "SELECT COUNT(*) FROM genotypes WHERE dosage IS NOT NULL AND rsid LIKE 'rs%'"
    ).fetchone()[0]
    print(f"[PRS-scorer] Allocating arrays for {n_total:,} genotypes...", flush=True)

    rsids = np.empty(n_total, dtype=np.int64)
    doses = np.empty(n_total, dtype=np.float64)  # float64 now!
    refs = np.empty(n_total, dtype=np.uint8)
    alts = np.empty(n_total, dtype=np.uint8)

    i = 0
    for rsid, dosage, ref, alt in db.execute(
        "SELECT rsid, dosage, ref, alt FROM genotypes WHERE dosage IS NOT NULL AND rsid LIKE 'rs%'"
    ):
        try:
            rsids[i] = int(rsid[2:])
            doses[i] = dosage
            refs[i] = BASE_ENCODE.get((ref or '').upper(), 4)
            alts[i] = BASE_ENCODE.get((alt or '').upper(), 4)
            i += 1
        except (ValueError, IndexError):
            pass
        if i % 5000000 == 0 and i > 0:
            print(f"[PRS-scorer]   {i:,}/{n_total:,} loaded...", flush=True)
    db.close()

    if i < n_total:
        rsids, doses, refs, alts = rsids[:i], doses[:i], refs[:i], alts[:i]

    print(f"[PRS-scorer] Sorting {i:,} entries...", flush=True)
    order = np.argsort(rsids)
    rsids, doses, refs, alts = rsids[order], doses[order], refs[order], alts[order]

    mem_mb = (rsids.nbytes + doses.nbytes + refs.nbytes + alts.nbytes) / 1e6
    print(f"[PRS-scorer] Loaded {i:,} genotypes into {mem_mb:.0f} MB ({time.time()-t0:.1f}s)", flush=True)
    return rsids, doses, refs, alts


def load_eur_freqs_numpy():
    """Load EUR allele frequencies into sorted numpy arrays for ambiguous SNP resolution."""
    t0 = time.time()
    if not os.path.exists(EUR_FREQ_FILE):
        print(f"[PRS-scorer] WARNING: No EUR freq file, ambiguous SNPs will be skipped", flush=True)
        return np.array([], dtype=np.int64), np.array([], dtype=np.uint8), np.array([], dtype=np.float32)

    print(f"[PRS-scorer] Loading EUR allele frequencies...", flush=True)
    with gzip.open(EUR_FREQ_FILE, "rt") as f:
        eur_data = json.load(f)

    n = len(eur_data)
    rsids = np.empty(n, dtype=np.int64)
    alleles = np.empty(n, dtype=np.uint8)
    freqs = np.empty(n, dtype=np.float32)

    i = 0
    for rsid_str, val in eur_data.items():
        if not rsid_str.startswith('rs'):
            continue
        try:
            rsids[i] = int(rsid_str[2:])
            alleles[i] = BASE_ENCODE.get(val['a'].upper(), 4)
            freqs[i] = val['f']
            i += 1
        except:
            pass

    rsids, alleles, freqs = rsids[:i], alleles[:i], freqs[:i]
    order = np.argsort(rsids)
    rsids, alleles, freqs = rsids[order], alleles[order], freqs[order]

    print(f"[PRS-scorer] {i:,} EUR freq variants ({time.time()-t0:.1f}s)", flush=True)
    return rsids, alleles, freqs


def parse_score_file(fpath):
    """Parse and score one PGS file with EUR freq-based ambiguous SNP resolution.

    For non-ambiguous SNPs (A/G, A/C, T/G, T/C): standard strand-flip logic.
    For ambiguous SNPs (A/T, C/G): use EUR allele frequency to determine correct
    allele orientation. Skip if EUR freq unavailable or freq ≈ 0.5 (unresolvable).

    Returns (score, matched, total, n_ambig_resolved, n_ambig_skipped) or None."""

    rsid_ints_list = []
    ea_codes_list = []
    weights_list = []

    with open(fpath) as f:
        for line in f:
            if line[0] != 'r':
                continue
            tab1 = line.index('\t')
            tab2 = line.index('\t', tab1 + 1)
            rsid_str = line[:tab1]
            ea_char = line[tab1+1:tab2].strip().upper()
            try:
                rsid_int = int(rsid_str[2:])
                tab3 = line.find('\t', tab2 + 1)
                w_str = line[tab2+1:tab3] if tab3 > 0 else line[tab2+1:].strip()
                weight = float(w_str)
            except:
                continue

            ea_code = BASE_ENCODE.get(ea_char, 4)
            if ea_code == 4:
                continue  # skip non-ACGT alleles
            rsid_ints_list.append(rsid_int)
            ea_codes_list.append(ea_code)
            weights_list.append(weight)

    total = len(rsid_ints_list)
    if total == 0:
        return None

    rsid_ints = np.array(rsid_ints_list, dtype=np.int64)
    ea_codes = np.array(ea_codes_list, dtype=np.uint8)
    weights = np.array(weights_list, dtype=np.float64)
    del rsid_ints_list, ea_codes_list, weights_list

    # Batch lookup in genotype arrays
    positions = np.searchsorted(_RSID_ARR, rsid_ints)
    positions_clamped = np.minimum(positions, len(_RSID_ARR) - 1)
    found = _RSID_ARR[positions_clamped] == rsid_ints

    if not np.any(found):
        return None

    # Extract matched data
    m_pos = positions_clamped[found]
    m_ea = ea_codes[found]
    m_w = weights[found]
    m_dos = _DOS_ARR[m_pos]  # already float64
    m_ref = _REF_ARR[m_pos]
    m_alt = _ALT_ARR[m_pos]
    m_rsid = rsid_ints[found]  # keep rsids for EUR freq lookup

    ea_comp = COMP_TABLE[m_ea]

    # Standard allele matching
    match_alt = (m_ea == m_alt) | (ea_comp == m_alt)
    match_ref = (m_ea == m_ref) | (ea_comp == m_ref)

    # Identify ambiguous SNPs (A/T or C/G pairs)
    is_ambig = (
        ((m_ref == 0) & (m_alt == 3)) | ((m_ref == 3) & (m_alt == 0)) |  # A/T
        ((m_ref == 1) & (m_alt == 2)) | ((m_ref == 2) & (m_alt == 1))    # C/G
    )

    # --- NON-AMBIGUOUS SNPs: standard logic ---
    non_ambig = ~is_ambig & (match_alt | match_ref)
    na_dos = m_dos[non_ambig].copy()
    na_flip = match_ref[non_ambig] & ~match_alt[non_ambig]
    na_dos[na_flip] = 2.0 - na_dos[na_flip]
    na_w = m_w[non_ambig]

    # --- AMBIGUOUS SNPs: EUR freq-based resolution ---
    n_ambig_resolved = 0
    n_ambig_skipped = 0
    ambig_dos_list = []
    ambig_w_list = []

    if np.any(is_ambig) and len(_EUR_RSID) > 0:
        ambig_idx = np.where(is_ambig)[0]
        ambig_rsids = m_rsid[ambig_idx]

        # Batch EUR freq lookup
        eur_pos = np.searchsorted(_EUR_RSID, ambig_rsids)
        eur_pos_clamped = np.minimum(eur_pos, len(_EUR_RSID) - 1)
        eur_found = _EUR_RSID[eur_pos_clamped] == ambig_rsids

        for j in range(len(ambig_idx)):
            i_orig = ambig_idx[j]
            if not eur_found[j]:
                n_ambig_skipped += 1
                continue

            eur_allele = _EUR_ALLELE[eur_pos_clamped[j]]
            eur_freq = float(_EUR_FREQ[eur_pos_clamped[j]])

            # Skip if MAF > 0.4 (freq too close to 0.5 to resolve)
            if 0.4 < eur_freq < 0.6:
                n_ambig_skipped += 1
                continue

            ea = m_ea[i_orig]
            ref = m_ref[i_orig]
            alt = m_alt[i_orig]
            dos = m_dos[i_orig]

            # Determine if ea corresponds to alt or ref using EUR freq
            # EUR says: eur_allele has frequency eur_freq
            # If eur_allele == alt → alt_freq = eur_freq
            #   If ea == alt → ea IS alt → no flip
            #   Else (ea == ref via complement) → ea IS ref → flip
            # If eur_allele == ref → ref_freq = eur_freq → alt_freq = 1-eur_freq
            #   If ea == ref → ea IS ref → flip
            #   Else (ea == alt via complement) → ea IS alt → no flip
            should_flip = False
            if eur_allele == alt:
                # EUR tracks alt allele
                if ea != alt:  # ea must be ref (complement of alt for ambig)
                    should_flip = True
            elif eur_allele == ref:
                # EUR tracks ref allele
                if ea == ref:  # ea directly matches ref
                    should_flip = True
                # else: ea is complement of ref = alt → no flip
            else:
                # EUR allele doesn't match either DB allele (rare)
                n_ambig_skipped += 1
                continue

            d = (2.0 - dos) if should_flip else dos
            ambig_dos_list.append(d)
            ambig_w_list.append(m_w[i_orig])
            n_ambig_resolved += 1
    else:
        n_ambig_skipped = int(np.sum(is_ambig))

    # Combine non-ambiguous + resolved ambiguous
    if ambig_dos_list:
        all_dos = np.concatenate([na_dos, np.array(ambig_dos_list)])
        all_w = np.concatenate([na_w, np.array(ambig_w_list)])
    else:
        all_dos = na_dos
        all_w = na_w

    matched = len(all_dos)
    if matched == 0:
        return None

    score = float(np.dot(all_dos, all_w))
    return (score, matched, total, n_ambig_resolved, n_ambig_skipped)


def normal_cdf(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def load_ld_inflation():
    ld_inflation = {}
    try:
        with open(METADATA_CSV, 'r') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if not row or len(row) < 9:
                    continue
                pgs_id = row[0].strip()
                if not pgs_id.startswith('PGS'):
                    continue
                line_lower = ','.join(row).lower()
                n_var = 0
                try:
                    n_var = int(row[8]) if row[8].strip() else 0
                except (ValueError, IndexError):
                    pass
                if re.search(r'ldpred|prs-cs|prscs|bayesian|lassosum|sblup|sbayesr|sbayess|ldak|megaprs', line_lower):
                    ld_inflation[pgs_id] = 1.0
                elif re.search(r'clumping|c\+t|p\+t|pruning|thresholding', line_lower):
                    ld_inflation[pgs_id] = 1.15
                else:
                    if n_var > 100000:
                        ld_inflation[pgs_id] = 1.4
                    elif re.search(r'genome-wide|genomewide', line_lower):
                        ld_inflation[pgs_id] = 1.3
                    else:
                        ld_inflation[pgs_id] = 1.25
    except Exception as e:
        print(f"[PRS-scorer] Warning: LD inflation load failed: {e}", flush=True)
    return ld_inflation


def calibrate_scores(results, dists):
    """Calibrate using pre-computed EUR freq stats (preferred) or distribution fallback."""
    if not results:
        return results

    ld_inflation = load_ld_inflation()
    freq_stats = {}
    if os.path.exists(FREQ_STATS_CACHE):
        with open(FREQ_STATS_CACHE) as f:
            freq_stats = json.load(f)
        print(f"[PRS-scorer] Loaded freq stats cache ({len(freq_stats)} models)", flush=True)

    freq_scored = 0
    dist_scored = 0
    fallback = 0

    for s in results:
        pgs_id = s["pgs_id"]
        raw_score = s["raw_score"]
        matched = s.get("matched", 0)

        fs = freq_stats.get(pgs_id)
        if fs and fs.get("freq_matched", 0) >= 10 and matched >= 10:
            f_mean = fs["freq_mean"]
            f_var = fs["freq_var"]
            f_matched = fs["freq_matched"]
            scale = matched / f_matched
            expected_mean = f_mean * scale
            ld_factor = ld_inflation.get(pgs_id, 1.25)
            expected_std = math.sqrt(f_var * ld_factor * scale) if f_var > 0 else 1e-10
            z = (raw_score - expected_mean) / expected_std
            z = max(-6.0, min(6.0, z))
            pct = max(0.1, min(99.9, normal_cdf(z) * 100))
            s["z_score"] = round(z, 2)
            s["percentile"] = round(pct, 1)
            s["scoring_method"] = "freq_based"
            s["ld_factor"] = ld_factor
            freq_scored += 1
            continue

        dist = dists.get(pgs_id)
        if dist:
            ref_scores = dist.get("EUR", {}).get("scores", [])
            if ref_scores and len(ref_scores) >= 20:
                ref_mean = statistics.mean(ref_scores)
                ref_std = statistics.stdev(ref_scores) if len(ref_scores) > 1 else 1.0
                if ref_std == 0:
                    ref_std = 1e-10
                z = (raw_score - ref_mean) / ref_std
                ld_factor = ld_inflation.get(pgs_id, 1.25)
                z = z / math.sqrt(ld_factor)
                z = max(-6.0, min(6.0, z))
                pct = max(0.1, min(99.9, normal_cdf(z) * 100))
                s["z_score"] = round(z, 2)
                s["percentile"] = round(pct, 1)
                s["scoring_method"] = "dist_based"
                s["ld_factor"] = ld_factor
                dist_scored += 1
                continue
        fallback += 1

    print(f"[PRS-scorer] Calibrated: {freq_scored} freq-based, {dist_scored} dist-based, {fallback} uncalibrated", flush=True)
    return results


def load_efo_labels():
    """Load EFO trait labels from PGS metadata CSV. Returns dict: pgs_id -> efo_label."""
    efo_map = {}
    try:
        with open(METADATA_CSV, 'r') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if len(row) > 4:
                    pgs_id = row[0].strip()
                    efo_label = row[3].strip().lower()
                    if pgs_id.startswith('PGS') and efo_label:
                        efo_map[pgs_id] = efo_label
        print(f"[PRS-scorer] Loaded {len(efo_map)} EFO trait labels", flush=True)
    except Exception as e:
        print(f"[PRS-scorer] Warning: EFO labels load failed: {e}", flush=True)
    return efo_map


def build_ensemble(results):
    """Group calibrated scores by EFO trait, remove outliers, produce ensemble z-average.

    Uses PGS Catalog's standardized EFO trait labels for grouping (much better than
    string normalization — "type 2 diabetes mellitus" groups all 70+ T2D models).

    For traits with multiple models:
    1. Compute median z-score across models
    2. Remove models deviating > 1.0 from median (likely allele orientation errors)
    3. Weighted average of remaining z-scores (weight = coverage)
    4. Convert to percentile

    Returns list of ensemble scores."""
    efo_map = load_efo_labels()

    trait_groups = {}
    for s in results:
        if 'z_score' not in s:
            continue
        # Exclude models with <20 variants - not true polygenic scores
        if s.get('matched', s.get('total', 0)) < 20:
            continue
        pgs_id = s['pgs_id']
        # Use EFO label if available, otherwise fall back to trait name
        key = efo_map.get(pgs_id, '')
        if not key:
            key = s.get('trait', pgs_id).lower().strip()
        s['ensemble_trait_key'] = key
        trait_groups.setdefault(key, []).append(s)

    ensemble = []

    for key, models in trait_groups.items():
        if len(models) == 1:
            s = models[0].copy()
            s['ensemble'] = False
            s['n_models'] = 1
            ensemble.append(s)
            continue

        z_scores = np.array([m['z_score'] for m in models])
        coverages = np.array([m.get("coverage", 0.5) for m in models])

        median_z = float(np.median(z_scores))

        # Remove outliers: models deviating > 1.0 from median
        keep = np.abs(z_scores - median_z) <= 1.0
        if not np.any(keep):
            keep = np.ones(len(models), dtype=bool)  # fallback: keep all

        kept_z = z_scores[keep]
        kept_cov = coverages[keep]
        kept_models = [m for m, k in zip(models, keep) if k]

        # Equal-weight average (trimmed mean) — no single model dominates
        avg_z = float(np.mean(kept_z))
        avg_z = max(-6.0, min(6.0, avg_z))
        pct = max(0.1, min(99.9, normal_cdf(avg_z) * 100))

        # Use highest-coverage model's metadata
        best = max(kept_models, key=lambda m: m.get('coverage', 0))

        ens_entry = {
            'pgs_id': best['pgs_id'],
            'trait': best.get('trait', key),
            'category': best.get('category', 'other'),
            'raw_score': best.get('raw_score', 0),
            'z_score': round(avg_z, 2),
            'percentile': round(pct, 1),
            'percentile_corrected': round(pct, 1),
            'score_raw': best.get('raw_score', 0),
            'corrected_score': best.get('raw_score', 0),
            'ensemble': True,
            'n_models': int(np.sum(keep)),
            'n_models_total': len(models),
            'n_outliers_removed': int(np.sum(~keep)),
            'model_z_min': round(float(np.min(kept_z)), 2),
            'model_z_max': round(float(np.max(kept_z)), 2),
            'model_z_spread': round(float(np.max(kept_z) - np.min(kept_z)), 2),
            'best_pgs_id': best['pgs_id'],
            'all_pgs_ids': [m['pgs_id'] for m in kept_models],
            'scoring_method': 'ensemble',
            'ensemble_trait_key': key,
            'matched': best.get('matched', 0),
            'total': best.get('total', 0),
            'coverage': best.get('coverage', 0),
            'freq_matched': best.get('matched', 0),
            'ref_pop': 'EUR',
            'ref_n': best.get('ref_n', 4257),
        }
        ensemble.append(ens_entry)

    return ensemble


def main():
    global _RSID_ARR, _DOS_ARR, _REF_ARR, _ALT_ARR, _EUR_RSID, _EUR_ALLELE, _EUR_FREQ
    if len(sys.argv) < 3:
        print("Usage: prs-scorer.py <genotype-db> <output-json> [--name NAME]")
        sys.exit(1)

    geno_file = sys.argv[1]
    out_file = sys.argv[2]
    name = "unknown"
    if "--name" in sys.argv:
        ni = sys.argv.index("--name")
        if ni + 1 < len(sys.argv):
            name = sys.argv[ni + 1]

    if not geno_file.endswith('.db'):
        print(f"[PRS-scorer] ERROR: requires .db input", flush=True)
        sys.exit(1)

    checkpoint_file = out_file + ".checkpoint"
    t0 = time.time()

    # 1. Load genotypes (~600MB with float64)
    _RSID_ARR, _DOS_ARR, _REF_ARR, _ALT_ARR = load_genotypes_numpy(geno_file)

    # 2. Load EUR allele frequencies (~46MB) for ambiguous SNP resolution
    _EUR_RSID, _EUR_ALLELE, _EUR_FREQ = load_eur_freqs_numpy()

    # 3. Load distributions + metadata
    print(f"[PRS-scorer] Loading distributions...", flush=True)
    with gzip.open(DIST_FILE, "rt") as f:
        dists = json.load(f)
    trait_info = {}
    if os.path.exists(TRAIT_FILE):
        with open(TRAIT_FILE) as f:
            trait_info = json.load(f)
    pgs_meta = {}
    meta_path = "/opt/helix/data/pgs-metadata-lookup.json"
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            pgs_meta = json.load(f)

    # 4. Score files
    # Filter out genome-wide models (>100K variants) - imputation bias inflates their scores
    all_score_files = sorted([
        fn for fn in os.listdir(SCORE_DIR)
        if fn.endswith('.txt') and fn.replace('.txt', '') in dists
    ])
    score_files = []
    skipped_big = 0
    for fn in all_score_files:
        fpath_check = os.path.join(SCORE_DIR, fn)
        # Quick line count via file size heuristic: avg ~30 bytes/line, >100K lines ~ >3MB
        fsize = os.path.getsize(fpath_check)
        if fsize > 3_000_000:  # >3MB = likely >100K variants
            # Verify with actual line count (fast - just count newlines)
            with open(fpath_check, 'rb') as fcheck:
                n_lines = fcheck.read().count(b'\n')
            if n_lines > 100_000:
                skipped_big += 1
                continue
        score_files.append(fn)
    n_models = len(score_files)
    if skipped_big > 0:
        print(f"[PRS-scorer] Skipped {skipped_big} genome-wide models (>100K variants)", flush=True)

    # 5. Checkpoint
    all_results = []
    scored_ids = set()
    resume_from = 0
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file) as f:
                ckpt = json.load(f)
            all_results = ckpt.get("scores", [])
            scored_ids = set(r["pgs_id"] for r in all_results)
            resume_from = ckpt.get("next_batch", 0)
            print(f"[PRS-scorer] RESUMING: {len(all_results)} done, batch {resume_from}", flush=True)
        except:
            pass

    n_workers = min(32, (os.cpu_count() or 4) * 4)
    total_ambig_resolved = 0
    total_ambig_skipped = 0
    print(f"[PRS-scorer] {n_models} models, {len(_RSID_ARR):,} genotypes, "
          f"{len(_EUR_RSID):,} EUR freqs, {n_workers} processes", flush=True)

    t1 = time.time()
    for batch_start in range(resume_from, n_models, BATCH_SIZE):
        batch_files = score_files[batch_start:batch_start + BATCH_SIZE]

        with ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=multiprocessing.get_context('fork')
        ) as pool:
            futures = {}
            for fn in batch_files:
                pgs_id = fn.replace('.txt', '')
                if pgs_id in scored_ids:
                    continue
                fpath = os.path.join(SCORE_DIR, fn)
                futures[pool.submit(parse_score_file, fpath)] = pgs_id

            for fut in as_completed(futures):
                pgs_id = futures[fut]
                try:
                    result = fut.result()
                    if result is None:
                        continue
                    score, matched, total, n_ar, n_as = result
                    total_ambig_resolved += n_ar
                    total_ambig_skipped += n_as

                    dist = dists.get(pgs_id)
                    if not dist:
                        continue
                    ref_scores = dist.get("EUR", {}).get("scores", [])
                    if not ref_scores or len(ref_scores) < 20:
                        continue

                    info = trait_info.get(pgs_id, pgs_meta.get(pgs_id, {}))
                    coverage = matched / max(total, 1)

                    all_results.append({
                        "pgs_id": pgs_id,
                        "trait": info.get("trait", ""),
                        "category": info.get("category", "other"),
                        "raw_score": round(score, 6),
                        "z_score": 0,
                        "percentile": 50.0,
                        "matched": matched,
                        "total": total,
                        "coverage": round(coverage, 4),
                        "ref_pop": "EUR",
                        "ref_n": len(ref_scores),
                        "method": "numpy_v6",
                        "ambig_resolved": n_ar,
                        "ambig_skipped": n_as,
                    })
                    scored_ids.add(pgs_id)
                except Exception as e:
                    print(f"[PRS-scorer] Error {pgs_id}: {e}", flush=True)

        done = min(batch_start + BATCH_SIZE, n_models)
        elapsed = time.time() - t1
        rate = done / elapsed if elapsed > 0 else 0
        eta = (n_models - done) / rate if rate > 0 else 0
        print(f"[PRS-scorer] {done}/{n_models} ({rate:.1f}/s, {len(all_results)} scored, ETA {eta:.0f}s)", flush=True)

        with open(checkpoint_file, 'w') as f:
            json.dump({"scores": all_results, "next_batch": batch_start + BATCH_SIZE,
                        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ"), "name": name}, f)

    print(f"[PRS-scorer] Ambiguous SNPs: {total_ambig_resolved:,} resolved, {total_ambig_skipped:,} skipped", flush=True)

    # Calibrate individual scores
    all_results = calibrate_scores(all_results, dists)

    # Add server.js compatibility fields
    for s in all_results:
        s["percentile_corrected"] = s["percentile"]
        s["score_raw"] = s["raw_score"]
        s["corrected_score"] = s["raw_score"]
        s["freq_matched"] = s.get("matched", 0)

    # Build ensemble scores
    ensemble_results = build_ensemble(all_results)
    n_ensemble = sum(1 for s in ensemble_results if s.get('ensemble'))
    n_single = sum(1 for s in ensemble_results if not s.get('ensemble'))
    total_outliers = sum(s.get('n_outliers_removed', 0) for s in ensemble_results)
    print(f"[PRS-scorer] Ensemble: {n_ensemble} multi-model traits, {n_single} single-model, "
          f"{total_outliers} outlier models removed", flush=True)

    elevated = sum(1 for s in ensemble_results if s.get("percentile", 50) >= 80)
    protective = sum(1 for s in ensemble_results if s.get("percentile", 50) <= 20)
    print(f"[PRS-scorer]   {elevated} elevated ({100*elevated/len(ensemble_results):.1f}%), "
          f"{protective} protective ({100*protective/len(ensemble_results):.1f}%)", flush=True)

    # Sort by trait name
    all_results.sort(key=lambda s: s.get("trait", ""))
    ensemble_results.sort(key=lambda s: s.get("trait", ""))

    with open(out_file, 'w') as f:
        json.dump({
            "scores": all_results,
            "ensemble_scores": ensemble_results,
            "method": "numpy_v6",
            "ref_panel": "opensnp_hm3",
            "n_ref": 4257,
            "n_models": len(all_results),
            "n_ensemble_traits": len(ensemble_results),
            "ambig_resolved": total_ambig_resolved,
            "ambig_skipped": total_ambig_skipped,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        }, f)

    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    print(f"[PRS-scorer] DONE: {len(all_results)} individual + {len(ensemble_results)} ensemble "
          f"in {time.time()-t0:.0f}s total", flush=True)


if __name__ == "__main__":
    main()
