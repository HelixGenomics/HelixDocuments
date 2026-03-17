#!/usr/bin/env python3
"""
PGP PRS GPU Scoring Pipeline — Final Version

Architecture:
  Phase 0: Build rsid index (background thread, one-time)
  Phase 1: Collect valid rsIDs (64 workers querying by pgs_id, ~2 min)
  Phase 2: Parse VCFs with rsID filtering (16 workers, ~15 min)
  Phase 3: Batched GPU scoring (load batch → sparse matmul → next batch)
  Phase 4: Save scores + calibrate

Memory budget: ~250GB for 174 VCF genotype dicts (filtered to ~10M entries each)
GPU budget: ~6GB VRAM per batch (sparse weights + dense dosage matrix)

Build machine: 256 cores, 503GB RAM, RTX 3060 Ti 8GB
"""
import os, sys, json, gzip, sqlite3, time, gc, math, argparse, threading
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
import scipy.sparse as sp

try:
    import torch
    HAS_CUDA = torch.cuda.is_available()
    if HAS_CUDA:
        GPU_NAME = torch.cuda.get_device_name(0)
        GPU_MEM = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        print("GPU: %s (%.1f GB)" % (GPU_NAME, GPU_MEM), flush=True)
    else:
        print("No CUDA GPU, using CPU sparse matmul", flush=True)
except ImportError:
    HAS_CUDA = False
    print("PyTorch not installed, using CPU sparse matmul", flush=True)

# ═══════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════

MEGA_DB = "/workspace/pgp-pipeline/pgs-mega.db"
VCF_DIR = "/workspace/pgp-pipeline/imputed-vcfs"
QC_FILE = "/workspace/pgp-pipeline/vcf-qc-results.json"
TRAIT_MAP_FILE = "/workspace/pgp-pipeline/pgs-trait-map.json"
OUTPUT_FILE = "/workspace/pgp-pipeline/pgp-prs-scores-gpu.json.gz"
DISTS_FILE = "/workspace/pgp-pipeline/opensnp-dists.json"
PROD_URL = "http://47.165.156.55:44255"
API_KEY = "hx-VhvdWnvq_dTCJI8r6CYRmzg5zY_EVY1Z"

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Module-level globals for fork-based COW sharing
_VALID_RSIDS = set()
_MEGA_DB_PATH = MEGA_DB


def post_progress(msg):
    try:
        import urllib.request
        data = json.dumps({"status": msg, "phase": "gpu-scoring"}).encode()
        req = urllib.request.Request(
            "%s/api/pipeline-status" % PROD_URL, data=data,
            headers={"Content-Type": "application/json", "X-Api-Key": API_KEY},
            method="POST")
        urllib.request.urlopen(req, timeout=3)
    except:
        pass


# ═══════════════════════════════════════════
# Phase 0: Build rsid index (background)
# ═══════════════════════════════════════════

def build_rsid_index():
    """Create rsid index on mega DB. Runs in background thread."""
    db = sqlite3.connect(MEGA_DB)
    # Check if index already exists
    existing = db.execute(
        "SELECT name FROM sqlite_master WHERE type='index' AND name='idx_rsid'"
    ).fetchone()
    if existing:
        print("  [IDX] rsid index already exists, skipping", flush=True)
        db.close()
        return

    print("  [IDX] Creating rsid index (background, may take 30-60 min)...", flush=True)
    db.execute("PRAGMA threads=64")
    db.execute("PRAGMA cache_size=-4000000")  # 4GB cache
    t0 = time.time()
    try:
        db.execute("CREATE INDEX idx_rsid ON pgs_variants(rsid)")
        db.commit()
        print("  [IDX] rsid index created in %.0fs" % (time.time() - t0), flush=True)
    except Exception as e:
        print("  [IDX] Index creation failed: %s" % e, flush=True)
    db.close()


# ═══════════════════════════════════════════
# Phase 1: Collect valid rsIDs (parallel)
# ═══════════════════════════════════════════

def _collect_rsids_worker(pgs_id):
    """Worker: collect distinct rsIDs for one PGS model. Uses pgs_id index = fast."""
    try:
        db = sqlite3.connect("file://%s?mode=ro" % _MEGA_DB_PATH, uri=True)
        db.execute("PRAGMA mmap_size=268435456")  # 256MB mmap per worker
        rsids = set()
        for row in db.execute(
            "SELECT rsid FROM pgs_variants WHERE pgs_id=?", (pgs_id,)
        ):
            if row[0]:
                rsids.add(row[0])
        db.close()
        return pgs_id, rsids
    except Exception as e:
        return pgs_id, set()


RSID_FILE = os.path.join(os.path.dirname(MEGA_DB), "valid-rsids.txt")


def collect_valid_rsids(n_workers=64):
    """Collect all unique rsIDs. Prefers cached file, falls back to DB query."""
    global _VALID_RSIDS
    t0 = time.time()

    # Fast path: load from pre-dumped file
    if os.path.exists(RSID_FILE):
        print("  Loading rsIDs from %s..." % RSID_FILE, flush=True)
        rsids = set()
        with open(RSID_FILE) as f:
            for line in f:
                rs = line.strip()
                if rs:
                    rsids.add(rs)
        _VALID_RSIDS = rsids
        print("  Loaded %d unique rsIDs from file in %.0fs" % (len(rsids), time.time() - t0), flush=True)
        return rsids

    # DB path: use rsid index with large mmap
    print("  No rsID cache file, querying DB...", flush=True)
    db = sqlite3.connect("file://%s?mode=ro" % MEGA_DB, uri=True)
    db.execute("PRAGMA mmap_size=200000000000")  # 200GB
    db.execute("PRAGMA cache_size=-8000000")      # 8GB cache
    rsids = set()
    for row in db.execute("SELECT DISTINCT rsid FROM pgs_variants"):
        if row[0]:
            rsids.add(row[0])
        if len(rsids) % 2000000 == 0 and len(rsids) > 0:
            print("    %d unique rsIDs, %.0fs..." % (len(rsids), time.time() - t0), flush=True)
    db.close()

    # Save for next time
    print("  Saving rsIDs to %s..." % RSID_FILE, flush=True)
    with open(RSID_FILE, 'w') as f:
        for rs in rsids:
            f.write(rs + '\n')

    _VALID_RSIDS = rsids
    print("  Collected %d unique rsIDs in %.0fs" % (len(rsids), time.time() - t0), flush=True)
    return rsids


# ═══════════════════════════════════════════
# Phase 2: Parse VCFs (CPU parallel)
# ═══════════════════════════════════════════

def _parse_vcf_worker(task):
    """Parse VCF → {rsid: (dosage, ref, alt)}. Only stores rsIDs in _VALID_RSIDS."""
    sid, vcf_path = task
    t0 = time.time()
    genotypes = {}
    n_vars = 0
    n_matched = 0
    valid = _VALID_RSIDS  # local ref for speed

    try:
        opener = gzip.open if vcf_path.endswith('.gz') else open
        with opener(vcf_path, 'rt') as f:
            for line in f:
                if line[0] == '#':
                    continue
                n_vars += 1
                parts = line.split('\t', 10)
                if len(parts) < 10:
                    continue

                rsid = parts[2]

                # Skip if not in mega DB (saves ~2/3 memory)
                if rsid not in valid:
                    continue

                ref = sys.intern(parts[3].upper())
                alt_raw = parts[4].upper()
                # Handle multi-allelic: take first alt
                alt = sys.intern(alt_raw.split(',')[0] if ',' in alt_raw else alt_raw)

                gt_field = parts[9]
                # Extract GT portion (before first :)
                gt_only = gt_field.split(':', 1)[0]
                sep = '/' if '/' in gt_only else '|'
                gt_parts = gt_only.split(sep)
                try:
                    a1 = int(gt_parts[0])
                    a2 = int(gt_parts[1]) if len(gt_parts) > 1 else a1
                except ValueError:
                    continue

                # DS tag for imputed dosages (more accurate than hard GT)
                dosage = None
                if ':' in gt_field and 'DS' in parts[8]:
                    fmt_fields = parts[8].split(':')
                    gt_values = gt_field.split(':')
                    for fi, fn in enumerate(fmt_fields):
                        if fn == 'DS' and fi < len(gt_values):
                            try:
                                dosage = float(gt_values[fi])
                            except ValueError:
                                pass
                            break

                if dosage is None:
                    dosage = float(a1 + a2)

                n_matched += 1
                genotypes[rsid] = (dosage, ref, alt)

    except Exception as e:
        print("  ERROR parsing %s: %s" % (vcf_path, e), flush=True)

    return sid, genotypes, n_vars, n_matched, time.time() - t0


def get_vcf_files():
    """Get QC-passed VCF files."""
    qc_passed = set()
    if os.path.exists(QC_FILE):
        qc = json.load(open(QC_FILE))
        if isinstance(qc, list):
            for entry in qc:
                if entry.get('status') == 'PASS':
                    qc_passed.add(entry.get('sample_id', ''))
        elif isinstance(qc, dict):
            for sid, info in qc.items():
                if info.get('passed', False) or info.get('status') == 'PASS':
                    qc_passed.add(sid)

    vcf_files = {}
    for f in sorted(os.listdir(VCF_DIR)):
        if f.endswith('-imputed.vcf.gz'):
            sid = f.replace('-imputed.vcf.gz', '')
            if qc_passed and sid not in qc_passed:
                continue
            vcf_files[sid] = os.path.join(VCF_DIR, f)
        elif f.endswith('.vcf.gz'):
            sid = f.replace('.vcf.gz', '')
            if qc_passed and sid not in qc_passed:
                continue
            vcf_files[sid] = os.path.join(VCF_DIR, f)

    return vcf_files


# ═══════════════════════════════════════════
# Phase 3: Batched GPU scoring
# ═══════════════════════════════════════════

def _load_model_batch_worker(pgs_id):
    """Worker: load one model's variants from mega DB."""
    try:
        db = sqlite3.connect("file://%s?mode=ro" % _MEGA_DB_PATH, uri=True)
        db.execute("PRAGMA mmap_size=268435456")
        results = []
        for rsid, ea, oa, w in db.execute(
            "SELECT rsid, effect_allele, other_allele, effect_weight "
            "FROM pgs_variants WHERE pgs_id=?", (pgs_id,)
        ):
            if rsid and w is not None:
                try:
                    wf = float(w)
                    if abs(wf) > 1e-15:
                        results.append((pgs_id, rsid, ea or "", oa or "", wf))
                except (ValueError, TypeError):
                    pass
        db.close()
        return results
    except Exception as e:
        return []


def create_batches(model_sizes, max_millions=150):
    """Split models into batches with max total variant count per batch."""
    sorted_models = sorted(model_sizes.items(), key=lambda x: -x[1])
    batches = []
    current_batch = []
    current_total = 0
    max_total = max_millions * 1_000_000

    for pgs_id, n_vars in sorted_models:
        if current_total + n_vars > max_total and current_batch:
            batches.append(current_batch)
            current_batch = []
            current_total = 0
        current_batch.append(pgs_id)
        current_total += n_vars

    if current_batch:
        batches.append(current_batch)

    return batches


def get_model_sizes():
    """Get variant counts per model from pgs_meta."""
    db = sqlite3.connect("file://%s?mode=ro" % MEGA_DB, uri=True)
    sizes = {}
    for pid, cnt in db.execute("SELECT pgs_id, variant_count FROM pgs_meta"):
        sizes[pid] = cnt
    db.close()
    return sizes


def match_dosage(vcf_dosage, vcf_ref, vcf_alt, model_ea, model_oa):
    """Match alleles between VCF and model. Returns adjusted dosage or None."""
    # Direct match: effect allele = alt
    if model_ea == vcf_alt:
        return vcf_dosage
    # Flipped: effect allele = ref
    if model_ea == vcf_ref:
        return 2.0 - vcf_dosage
    # Complement strand
    comp_alt = COMPLEMENT.get(vcf_alt, '')
    comp_ref = COMPLEMENT.get(vcf_ref, '')
    if model_ea == comp_alt:
        return vcf_dosage
    if model_ea == comp_ref:
        return 2.0 - vcf_dosage
    return None


def score_batch(batch_pgs_ids, sample_genotypes, n_load_workers, use_gpu, gpu_budget_gb):
    """
    Load one batch of models, build sparse matrix, score all samples.
    Returns: {sid: {pgs_id: raw_score}}
    """
    n_models = len(batch_pgs_ids)
    t0 = time.time()

    # Step 1: Load all variants for this batch (parallel by model)
    print("    Loading %d models..." % n_models, flush=True)
    all_variants = []  # list of (pgs_id, rsid, ea, oa, weight)

    with Pool(n_load_workers) as pool:
        for chunk in pool.imap_unordered(_load_model_batch_worker, batch_pgs_ids, chunksize=2):
            all_variants.extend(chunk)

    load_time = time.time() - t0
    print("    Loaded %d variant entries in %.0fs" % (len(all_variants), load_time), flush=True)

    if not all_variants:
        return {}

    # Step 2: Build sparse weight matrix (n_models × n_unique_rsids)
    t1 = time.time()
    rsid_to_col = {}
    allele_map = {}  # col -> (ea, oa)
    model_to_row = {pid: i for i, pid in enumerate(batch_pgs_ids)}
    col_counter = 0

    rows = []
    cols = []
    vals = []

    for pgs_id, rsid, ea, oa, weight in all_variants:
        if rsid not in rsid_to_col:
            rsid_to_col[rsid] = col_counter
            allele_map[col_counter] = (ea, oa)
            col_counter += 1

        row = model_to_row.get(pgs_id)
        if row is not None:
            rows.append(row)
            cols.append(rsid_to_col[rsid])
            vals.append(weight)

    n_cols = col_counter
    del all_variants
    # NOTE: No gc.collect() — heap is 200GB+, traversal takes 5-15 min per call

    # Build scipy sparse CSR
    weight_csr = sp.csr_matrix(
        (np.array(vals, dtype=np.float32),
         (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32))),
        shape=(n_models, n_cols)
    )
    del rows, cols, vals

    build_time = time.time() - t1
    print("    Matrix: %d models × %d rsIDs, %d nnz (%.0fs)" % (
        n_models, n_cols, weight_csr.nnz, build_time), flush=True)

    # Step 3: Build dosage matrix and score
    t2 = time.time()
    sample_ids = sorted(sample_genotypes.keys())
    n_samples = len(sample_ids)

    # Estimate memory for dense dosage matrix
    dense_gb = n_cols * n_samples * 4 / (1024**3)
    sparse_gb = weight_csr.nnz * 12 / (1024**3)  # val + col + rowptr

    # Decide: GPU batch matmul or CPU per-sample
    if use_gpu and (dense_gb + sparse_gb) < gpu_budget_gb:
        # ── GPU path: batch all samples in one matmul ──
        print("    GPU batch matmul (%.1f GB dense + %.1f GB sparse)..." % (
            dense_gb, sparse_gb), flush=True)

        dosage_matrix = np.zeros((n_cols, n_samples), dtype=np.float32)

        for si, sid in enumerate(sample_ids):
            genos = sample_genotypes[sid]
            for rsid, col in rsid_to_col.items():
                g = genos.get(rsid)
                if g is None:
                    continue
                dosage, ref, alt = g
                ea, oa = allele_map[col]
                matched = match_dosage(dosage, ref, alt, ea, oa)
                if matched is not None:
                    dosage_matrix[col, si] = matched

        # Transfer to GPU and compute
        W_torch = torch.sparse_csr_tensor(
            torch.tensor(weight_csr.indptr, dtype=torch.int64),
            torch.tensor(weight_csr.indices, dtype=torch.int64),
            torch.tensor(weight_csr.data, dtype=torch.float32),
            size=(n_models, n_cols)
        ).cuda()

        D_torch = torch.tensor(dosage_matrix, dtype=torch.float32).cuda()
        del dosage_matrix

        scores_gpu = torch.sparse.mm(W_torch, D_torch)  # (n_models × n_samples)
        scores_np = scores_gpu.cpu().numpy()

        del W_torch, D_torch, scores_gpu
        torch.cuda.empty_cache()

    else:
        # ── CPU path: scipy sparse matmul, chunked if needed ──
        if dense_gb > 30:
            # Too large for one shot, process samples in chunks
            print("    CPU chunked scoring (%d samples, %.1f GB per chunk)..." % (
                n_samples, n_cols * 4 / (1024**3)), flush=True)
            scores_np = np.zeros((n_models, n_samples), dtype=np.float32)
            chunk_size = max(1, int(30 * (1024**3) / (n_cols * 4)))  # ~30GB chunks

            for si_start in range(0, n_samples, chunk_size):
                si_end = min(si_start + chunk_size, n_samples)
                chunk_sids = sample_ids[si_start:si_end]
                dosage_chunk = np.zeros((n_cols, len(chunk_sids)), dtype=np.float32)

                for ci, sid in enumerate(chunk_sids):
                    genos = sample_genotypes[sid]
                    for rsid, col in rsid_to_col.items():
                        g = genos.get(rsid)
                        if g is None:
                            continue
                        dosage, ref, alt = g
                        ea, oa = allele_map[col]
                        matched = match_dosage(dosage, ref, alt, ea, oa)
                        if matched is not None:
                            dosage_chunk[col, ci] = matched

                scores_np[:, si_start:si_end] = (weight_csr @ dosage_chunk).astype(np.float32)
                del dosage_chunk
        else:
            print("    CPU batch scoring (%.1f GB dense matrix)..." % dense_gb, flush=True)
            dosage_matrix = np.zeros((n_cols, n_samples), dtype=np.float32)

            for si, sid in enumerate(sample_ids):
                genos = sample_genotypes[sid]
                for rsid, col in rsid_to_col.items():
                    g = genos.get(rsid)
                    if g is None:
                        continue
                    dosage, ref, alt = g
                    ea, oa = allele_map[col]
                    matched = match_dosage(dosage, ref, alt, ea, oa)
                    if matched is not None:
                        dosage_matrix[col, si] = matched

            scores_np = (weight_csr @ dosage_matrix).astype(np.float32)
            del dosage_matrix

    score_time = time.time() - t2

    # Step 4: Package results
    batch_results = {}
    for si, sid in enumerate(sample_ids):
        batch_results[sid] = {}
        for mi, pgs_id in enumerate(batch_pgs_ids):
            score = float(scores_np[mi, si])
            if abs(score) > 1e-15:
                batch_results[sid][pgs_id] = score

    del weight_csr, scores_np, rsid_to_col, allele_map

    total_time = time.time() - t0
    n_scored = sum(len(v) for v in batch_results.values())
    print("    Scored: %d non-zero entries | Load %.0fs + Build %.0fs + Score %.0fs = %.0fs total" % (
        n_scored, load_time, build_time, score_time, total_time), flush=True)

    return batch_results


# ═══════════════════════════════════════════
# Phase 4: Calibration
# ═══════════════════════════════════════════

def norm_cdf(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def calibrate_scores(raw_scores, trait_map):
    """Convert raw scores to percentiles using OpenSNP distributions."""
    dists = {}
    if os.path.exists(DISTS_FILE):
        with open(DISTS_FILE) as f:
            dists_data = json.load(f)
        dists = dists_data.get("distributions", dists_data)
        print("  Loaded distributions for %d models" % len(dists), flush=True)

    calibrated = {}
    for sid, scores in raw_scores.items():
        calibrated[sid] = {}
        for pgs_id, raw_score in scores.items():
            meta = trait_map.get(pgs_id, {})
            trait = meta.get("trait", pgs_id)
            category = meta.get("category", "other")

            entry = {
                "score": raw_score,
                "trait": trait,
                "category": category,
            }

            # OpenSNP calibration
            dist = dists.get(pgs_id, {})
            eur = dist.get("EUR", dist.get("ALL", {}))
            if eur and "mean" in eur and "std" in eur:
                mean = eur["mean"]
                std = eur["std"]
                if std > 1e-15:
                    z = (raw_score - mean) / std
                    pctl = norm_cdf(z) * 100
                    pctl = max(0.1, min(99.9, pctl))
                    entry["percentile"] = round(pctl, 2)
                    entry["z_score"] = round(z, 3)

            calibrated[sid][pgs_id] = entry

    return calibrated


# ═══════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════

def main():
    global _VALID_RSIDS
    parser = argparse.ArgumentParser(description="PGP PRS GPU Scoring Pipeline")
    parser.add_argument("--vcf-workers", type=int, default=16,
                        help="CPU workers for VCF parsing (default: 16)")
    parser.add_argument("--load-workers", type=int, default=64,
                        help="Workers for model loading from DB (default: 64)")
    parser.add_argument("--batch-size", type=int, default=150,
                        help="Max variants per batch in millions (default: 150)")
    parser.add_argument("--gpu-budget", type=float, default=6.0,
                        help="GPU VRAM budget in GB (default: 6.0)")
    parser.add_argument("--no-gpu", action="store_true",
                        help="Force CPU-only scoring")
    parser.add_argument("--build-index", action="store_true",
                        help="Build rsid index before scoring (foreground, ~30 min)")
    args = parser.parse_args()

    use_gpu = HAS_CUDA and not args.no_gpu

    print("=" * 60)
    print("PGP PRS SCORING PIPELINE")
    print("VCF workers: %d | Load workers: %d | Batch: %dM" % (
        args.vcf_workers, args.load_workers, args.batch_size))
    print("GPU: %s | Time: %s" % (
        'CUDA' if use_gpu else 'CPU',
        time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())))
    print("=" * 60, flush=True)
    t_start = time.time()

    # Load trait metadata
    trait_map = {}
    if os.path.exists(TRAIT_MAP_FILE):
        with open(TRAIT_MAP_FILE) as f:
            trait_map = json.load(f)
    print("Trait map: %d models" % len(trait_map), flush=True)

    # ── Phase 0: rsid index (background or foreground) ──
    if args.build_index:
        print("\n" + "=" * 60)
        print("PHASE 0: Building rsid index (foreground)")
        print("=" * 60, flush=True)
        build_rsid_index()
    else:
        # Start index creation in background (non-blocking)
        idx_thread = threading.Thread(target=build_rsid_index, daemon=True)
        idx_thread.start()
        print("  [IDX] rsid index creation started in background", flush=True)

    # ── Phase 1: Collect valid rsIDs ──
    print("\n" + "=" * 60)
    print("PHASE 1: Collect valid rsIDs from mega DB")
    print("=" * 60, flush=True)
    collect_valid_rsids(n_workers=args.load_workers)
    post_progress("Collected %d valid rsIDs" % len(_VALID_RSIDS))

    # ── Phase 2: Parse VCFs ──
    print("\n" + "=" * 60)
    print("PHASE 2: Parse VCFs (CPU parallel, rsID-filtered)")
    print("=" * 60, flush=True)

    vcf_files = get_vcf_files()
    print("  %d VCFs to parse" % len(vcf_files), flush=True)

    tasks = [(sid, path) for sid, path in sorted(vcf_files.items())]
    sample_genotypes = {}
    sample_stats = {}

    t_parse = time.time()
    done = 0
    with Pool(args.vcf_workers) as pool:
        for sid, genos, n_vars, n_matched, elapsed in pool.imap_unordered(
            _parse_vcf_worker, tasks, chunksize=1
        ):
            done += 1
            sample_genotypes[sid] = genos
            sample_stats[sid] = (n_vars, n_matched, len(genos))
            print("  [%d/%d] %s: %.1fM vars, %d matched (%.1f%%), %d stored, %.0fs" % (
                done, len(tasks), sid, n_vars / 1e6, n_matched,
                100 * n_matched / max(1, n_vars), len(genos), elapsed), flush=True)

            if done % 20 == 0:
                mem_gb = sum(sys.getsizeof(g) for g in sample_genotypes.values()) / (1024**3)
                post_progress("Parsed %d/%d VCFs" % (done, len(tasks)))

    parse_time = time.time() - t_parse
    total_genos = sum(len(g) for g in sample_genotypes.values())
    print("\n  Parsed %d VCFs in %.0fs (%dM total genotype entries)" % (
        len(sample_genotypes), parse_time, total_genos // 1_000_000), flush=True)
    post_progress("Parsed %d VCFs, %dM genotypes" % (len(sample_genotypes), total_genos // 1_000_000))

    # Free valid rsids set (no longer needed)
    _VALID_RSIDS = set()

    # Disable gc during scoring — traversing 233GB heap takes 15+ min per gc.collect()
    gc.disable()
    print("  gc disabled for scoring phase (heap too large for cyclic GC)", flush=True)

    # ── Phase 3: Batched GPU scoring ──
    print("\n" + "=" * 60)
    print("PHASE 3: Batched model loading + scoring")
    print("=" * 60, flush=True)

    model_sizes = get_model_sizes()
    batches = create_batches(model_sizes, max_millions=args.batch_size)
    print("  %d models → %d batches (max %dM variants each)" % (
        len(model_sizes), len(batches), args.batch_size), flush=True)

    # Show batch sizes
    for i, batch in enumerate(batches):
        batch_vars = sum(model_sizes.get(pid, 0) for pid in batch)
        print("    Batch %d: %d models, %dM variants" % (
            i + 1, len(batch), batch_vars // 1_000_000), flush=True)

    # Score each batch
    all_scores = defaultdict(dict)  # sid -> pgs_id -> raw_score
    t_score = time.time()

    for bi, batch_pgs_ids in enumerate(batches):
        batch_vars = sum(model_sizes.get(pid, 0) for pid in batch_pgs_ids)
        print("\n  ── Batch %d/%d: %d models, %dM variants ──" % (
            bi + 1, len(batches), len(batch_pgs_ids), batch_vars // 1_000_000), flush=True)

        batch_results = score_batch(
            batch_pgs_ids, sample_genotypes,
            n_load_workers=args.load_workers,
            use_gpu=use_gpu,
            gpu_budget_gb=args.gpu_budget
        )

        # Merge into all_scores
        for sid, scores in batch_results.items():
            all_scores[sid].update(scores)

        del batch_results
        if use_gpu:
            torch.cuda.empty_cache()

        post_progress("Scored batch %d/%d" % (bi + 1, len(batches)))

    score_time = time.time() - t_score
    total_scored = sum(len(v) for v in all_scores.values())
    print("\n  Scoring complete: %d samples × %d total score entries in %.0fs" % (
        len(all_scores), total_scored, score_time), flush=True)

    # Free genotype data and re-enable gc
    del sample_genotypes
    gc.enable()
    gc.collect()

    # ── Phase 4: Calibrate + Save ──
    print("\n" + "=" * 60)
    print("PHASE 4: Calibrate + Save")
    print("=" * 60, flush=True)

    # Convert to regular dict
    raw_scores = dict(all_scores)

    # Save raw scores
    raw_file = OUTPUT_FILE.replace('.json.gz', '-raw.json.gz')
    print("  Saving raw scores to %s..." % raw_file, flush=True)
    with gzip.open(raw_file, 'wt') as f:
        json.dump(raw_scores, f)
    raw_size = os.path.getsize(raw_file) / (1024**2)
    print("  Raw scores: %.1f MB" % raw_size, flush=True)

    # Calibrate
    calibrated = calibrate_scores(raw_scores, trait_map)

    # Add sample stats
    for sid in calibrated:
        if sid in sample_stats:
            n_vars, n_matched, n_stored = sample_stats.get(sid, (0, 0, 0))
            calibrated[sid]["_meta"] = {
                "vcf_variants": n_vars,
                "db_matched": n_matched,
                "genotypes_stored": n_stored,
                "models_scored": len(calibrated[sid]) - 1,
            }

    # Save calibrated
    print("  Saving calibrated scores to %s..." % OUTPUT_FILE, flush=True)
    with gzip.open(OUTPUT_FILE, 'wt') as f:
        json.dump(calibrated, f)
    cal_size = os.path.getsize(OUTPUT_FILE) / (1024**2)
    print("  Calibrated scores: %.1f MB" % cal_size, flush=True)

    # ── Summary ──
    total_time = time.time() - t_start
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print("  Samples: %d" % len(calibrated))
    print("  Models: %d" % len(model_sizes))
    print("  Total scores: %d" % total_scored)
    print("  Parse time: %.0fs" % parse_time)
    print("  Score time: %.0fs" % score_time)
    print("  Total time: %.0fs (%.1f min)" % (total_time, total_time / 60))
    print("  Output: %s" % OUTPUT_FILE, flush=True)

    post_progress("Pipeline complete: %d samples, %d models, %.0f min" % (
        len(calibrated), len(model_sizes), total_time / 60))

    # Quick validation: show score distribution for first sample
    first_sid = sorted(calibrated.keys())[0] if calibrated else None
    if first_sid:
        sample = calibrated[first_sid]
        pctls = [v.get("percentile") for v in sample.values()
                 if isinstance(v, dict) and v.get("percentile") is not None]
        if pctls:
            print("\n  Validation (%s): %d models with percentiles" % (first_sid, len(pctls)))
            print("  Mean: %.1f, Median: %.1f, Min: %.1f, Max: %.1f" % (
                sum(pctls) / len(pctls),
                sorted(pctls)[len(pctls) // 2],
                min(pctls), max(pctls)))


if __name__ == "__main__":
    main()
