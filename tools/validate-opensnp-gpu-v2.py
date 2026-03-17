#!/usr/bin/env python3
"""
OpenSNP Validation v4: GPU-accelerated scoring via PyTorch sparse matmul.
Fixes from v3:
  - BED lookup [2,255,1,0] counts A1 alleles (matches allele-aligned weight matrix v3)
  - IrisPlex SNP-based eye color prediction (6 key SNPs)
  - Weight matrix path updated to v3 (allele-aligned)
"""
import sqlite3, os, sys, time, json, gzip, re, csv
import numpy as np
from scipy import sparse, stats
from collections import defaultdict, Counter
import torch

BASE_DIR = "/workspace/pgp-pipeline"
MEGA_DB = f"{BASE_DIR}/pgs-mega.db"
OPENSNP_DIR = f"{BASE_DIR}/opensnp-frei2024"
HM3_PREFIX = f"{OPENSNP_DIR}/imputed/opensnp_hm3"
PHENO_FULL = f"{OPENSNP_DIR}/phenotypes_202402180100.csv"
PHENO_QC = f"{OPENSNP_DIR}/qc/pheno.csv"
WEIGHT_NPZ = f"{BASE_DIR}/opensnp-weight-matrix-v3.npz"
OUTPUT = f"{BASE_DIR}/opensnp-validation-results-v2.json"

GPU_BATCH = 200      # samples per GPU batch
WEIGHT_CHUNKS = 4    # split weight matrix columns into this many chunks

# IrisPlex 6 SNPs for eye color prediction
# Reference: Walsh et al. 2011, Forensic Science International: Genetics
IRISPLEX_SNPS = {
    "rs12913832": {"gene": "HERC2/OCA2", "blue_weight": 3.0},   # strongest predictor
    "rs1800407":  {"gene": "OCA2",        "blue_weight": 0.5},
    "rs12896399": {"gene": "SLC24A4",     "blue_weight": 0.3},
    "rs16891982": {"gene": "SLC45A2",     "blue_weight": 0.8},
    "rs1393350":  {"gene": "TYR",         "blue_weight": 0.2},
    "rs12203592": {"gene": "IRF4",        "blue_weight": 0.4},
}


def read_bed(bed_path, n_samples, n_snps):
    print(f"Reading BED: {n_snps:,} SNPs x {n_samples:,} samples...", flush=True)
    t0 = time.time()
    bytes_per_snp = (n_samples + 3) // 4
    with open(bed_path, "rb") as f:
        magic = f.read(3)
        assert magic[:2] == b"\x6c\x1b" and magic[2] == 1, "Invalid BED"
        raw = np.frombuffer(f.read(), dtype=np.uint8).reshape(n_snps, bytes_per_snp)

    # Count A1 alleles (PLINK default):
    #   00 = hom A1 -> 2 copies of A1
    #   01 = missing -> 255
    #   10 = het     -> 1 copy of A1
    #   11 = hom A2  -> 0 copies of A1
    lookup = np.array([2, 255, 1, 0], dtype=np.uint8)

    dose = np.empty((n_snps, bytes_per_snp * 4), dtype=np.uint8)
    dose[:, 0::4] = lookup[raw & 0x03]
    dose[:, 1::4] = lookup[(raw >> 2) & 0x03]
    dose[:, 2::4] = lookup[(raw >> 4) & 0x03]
    dose[:, 3::4] = lookup[(raw >> 6) & 0x03]
    dose = dose[:, :n_samples]
    missing = np.sum(dose == 255)
    dose[dose == 255] = 0
    print(f"  Decoded in {time.time()-t0:.1f}s, {missing:,} missing ({missing/(n_snps*n_samples)*100:.4f}%)", flush=True)
    return dose


def scipy_csc_to_torch_sparse(csc_mat, device='cuda'):
    """Convert scipy CSC sparse -> torch sparse CSR on GPU (transposed for matmul)."""
    csr = csc_mat.T.tocsr()
    crow = torch.tensor(csr.indptr.astype(np.int32), dtype=torch.int32, device=device)
    col = torch.tensor(csr.indices.astype(np.int32), dtype=torch.int32, device=device)
    vals = torch.tensor(csr.data.astype(np.float32), dtype=torch.float32, device=device)
    return torch.sparse_csr_tensor(crow, col, vals, size=csr.shape, device=device)


def gpu_score(dose, W_scipy, n_models, n_samples):
    """Score using GPU with chunked weight matrix."""
    print(f"\n  GPU: {torch.cuda.get_device_name(0)}, "
          f"VRAM: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB"
          if hasattr(torch.cuda.get_device_properties(0), 'total_mem')
          else f"\n  GPU: {torch.cuda.get_device_name(0)}, "
          f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB", flush=True)

    scores = np.zeros((n_samples, n_models), dtype=np.float32)
    chunk_size = (n_models + WEIGHT_CHUNKS - 1) // WEIGHT_CHUNKS

    for ci in range(WEIGHT_CHUNKS):
        col_start = ci * chunk_size
        col_end = min((ci + 1) * chunk_size, n_models)
        W_chunk = W_scipy[:, col_start:col_end]
        n_chunk_models = col_end - col_start

        print(f"  Chunk {ci+1}/{WEIGHT_CHUNKS}: models [{col_start}-{col_end}), "
              f"nnz={W_chunk.nnz:,}", flush=True)

        W_gpu = scipy_csc_to_torch_sparse(W_chunk, device='cuda')
        vram = torch.cuda.memory_allocated() / 1e9
        print(f"    Weight on GPU: {vram:.2f} GB", flush=True)

        for si in range(0, n_samples, GPU_BATCH):
            se = min(si + GPU_BATCH, n_samples)
            dose_batch = dose[:, si:se].astype(np.float32)
            dose_gpu = torch.tensor(dose_batch, device='cuda')
            result = torch.sparse.mm(W_gpu, dose_gpu)
            scores[si:se, col_start:col_end] = result.T.cpu().numpy()
            del dose_gpu, result

        del W_gpu
        torch.cuda.empty_cache()

    return scores


def parse_height(val):
    if not val or val == "-":
        return None
    val = val.lower().strip()
    if val in ["gg", "rs6060371", "ttcc", "n/a", "unknown", ""]:
        return None
    if val == ">200cm":
        return 205.0
    if val == "1m87":
        val = "187 cm"
    if "at max" in val:
        val = val.split("at max")[0].strip()
    cat = {
        "tall ( >180cm )": 185.0,
        "average ( 165cm < x < 180cm )": 172.5,
        "short (<165cm)": 160.0,
        "short (<165cm )": 160.0,
    }
    for k, v in cat.items():
        if val == k:
            return v
    m = re.search(r"(\d+\.?\d*)\s*cm", val)
    if m:
        cm = float(m.group(1))
        return cm if 120 < cm < 230 else None
    m = re.search(r"^(\d)[.,](\d+)\s*m\b", val)
    if m:
        cm = float(m.group(1)) * 100 + float(m.group(2).ljust(2, "0")[:2])
        return cm if 120 < cm < 230 else None
    m = re.search(r"(\d+)['\u2018\u2019]\s*(\d+\.?\d*)?[\"\u201d\u2033]?", val)
    if m:
        ft, inch = float(m.group(1)), float(m.group(2) or 0)
        if 3 <= ft <= 7 and inch < 12:
            cm = ft * 30.48 + inch * 2.54
            return cm if 120 < cm < 230 else None
    m = re.search(r"(\d+)\s*ft\.?\s*(\d+\.?\d*)?\s*(?:in)?", val)
    if m:
        ft, inch = float(m.group(1)), float(m.group(2) or 0)
        if 3 <= ft <= 7 and inch < 12:
            cm = ft * 30.48 + inch * 2.54
            return cm if 120 < cm < 230 else None
    m = re.match(r"^(\d+\.?\d*)$", val)
    if m:
        num = float(m.group(1))
        if 120 < num < 230:
            return num
        if 48 < num < 84:
            return num * 2.54
    return None


def clean_hair(x):
    x = x.lower().strip()
    if x == "-" or not x:
        return None
    if "blond" in x:
        return "blonde"
    if "brown" in x:
        return "brown"
    if "black" in x:
        return "black"
    if "red" in x or "auburn" in x or "ginger" in x:
        return "red"
    return "other"


def clean_eye(x):
    x = x.lower().strip()
    if x == "-" or not x:
        return None
    if "brown" in x and "green" not in x and "hazel" not in x:
        return "brown"
    if "blue" in x and "green" not in x:
        return "blue"
    if "green" in x and "brown" not in x and "blue" not in x and "hazel" not in x:
        return "green"
    if "hazel" in x:
        return "hazel"
    if "gray" in x or "grey" in x:
        return "blue"
    return "other"


def load_ancestry(pheno_qc_path, sample_ids):
    pc1_map = {}
    try:
        with open(pheno_qc_path) as f:
            header = f.readline().strip().split("\t")
            pc1_idx = header.index("PC1") if "PC1" in header else None
            for line in f:
                parts = line.strip().split("\t")
                sid = parts[0]
                if pc1_idx is not None and len(parts) > pc1_idx:
                    try:
                        pc1_map[sid] = float(parts[pc1_idx])
                    except ValueError:
                        pass
    except Exception as e:
        print(f"  Warning: could not load PCA from {pheno_qc_path}: {e}", flush=True)
        return {}

    ancestry = {}
    for sid in sample_ids:
        pc1 = pc1_map.get(sid)
        if pc1 is None:
            continue
        if pc1 > 0.04:
            ancestry[sid] = "EUR"
        elif pc1 < -0.01:
            ancestry[sid] = "non-EUR"
        else:
            ancestry[sid] = "admixed"
    return ancestry


def irisplex_eye_color(dose, bim_rsid_map, samples, phenotypes):
    """
    IrisPlex-like eye color prediction using 6 key SNPs.

    Since PGS Catalog has no eye color models, we use direct SNP-based prediction.
    rs12913832 (HERC2) is the strongest predictor of blue vs brown eyes.
    """
    print("\n" + "=" * 60)
    print("EYE COLOR VALIDATION (IrisPlex SNPs)")
    print("=" * 60, flush=True)

    # Find IrisPlex SNPs in BIM
    found_snps = {}
    for rsid, info in IRISPLEX_SNPS.items():
        if rsid in bim_rsid_map:
            found_snps[rsid] = {**info, "row": bim_rsid_map[rsid]}

    print(f"  IrisPlex SNPs found in BIM: {len(found_snps)}/{len(IRISPLEX_SNPS)}")
    for rsid, info in found_snps.items():
        print(f"    {rsid} ({info['gene']}): row {info['row']}")

    if not found_snps:
        print("  ERROR: No IrisPlex SNPs found in dataset!")
        return []

    # Get eye color phenotypes
    eye_idx_list, eye_labels = [], []
    for i, sid in enumerate(samples):
        if sid in phenotypes and "eye" in phenotypes[sid]:
            eye_idx_list.append(i)
            eye_labels.append(phenotypes[sid]["eye"])

    print(f"  Eye color samples: {len(eye_idx_list)}")
    print(f"  Distribution: {Counter(eye_labels)}")

    if len(eye_idx_list) < 30:
        print("  Too few samples for validation")
        return []

    # Build simple scoring: weighted sum of A1 dosage for IrisPlex SNPs
    # Higher score = more "blue eye" alleles
    irisplex_scores = np.zeros(len(eye_idx_list), dtype=np.float32)
    snp_details = {}

    for rsid, info in found_snps.items():
        row = info["row"]
        snp_dose = dose[row, eye_idx_list].astype(np.float32)

        # For rs12913832: A1 allele (usually G) is the blue-eye allele
        # We use the dosage directly (A1 count) weighted by importance
        irisplex_scores += snp_dose * info["blue_weight"]

        # Per-SNP stats
        mean_dose = np.mean(snp_dose)
        snp_details[rsid] = {
            "gene": info["gene"],
            "mean_dosage": round(float(mean_dose), 3),
            "weight": info["blue_weight"],
        }

    print(f"\n  IrisPlex composite score: mean={np.mean(irisplex_scores):.3f}, "
          f"std={np.std(irisplex_scores):.3f}")

    # Validate: compute AUC for each eye color
    from sklearn.metrics import roc_auc_score

    results = []
    for color in ["blue", "brown", "green", "hazel"]:
        binary = np.array([1 if e == color else 0 for e in eye_labels])
        n_pos = int(np.sum(binary))
        n_neg = len(binary) - n_pos
        if n_pos < 10 or n_neg < 10:
            print(f"  {color}: skipped (n_cases={n_pos})")
            continue

        try:
            # For blue: higher score should predict blue (positive AUC)
            # For brown: higher score should predict NOT brown (invert)
            if color in ["blue", "green"]:
                auc = roc_auc_score(binary, irisplex_scores)
            else:
                auc = roc_auc_score(binary, -irisplex_scores)

            # Also compute with max(auc, 1-auc) for direction-invariant metric
            auc_best = max(auc, 1.0 - auc)

            result = {
                "method": "irisplex",
                "target": color,
                "auc": round(float(auc_best), 4),
                "auc_directional": round(float(auc), 4),
                "n": len(eye_labels),
                "n_cases": n_pos,
                "n_snps_used": len(found_snps),
                "snps": list(found_snps.keys()),
            }
            results.append(result)
            print(f"  {color}: AUC={auc_best:.4f} (directional={auc:.4f}), "
                  f"cases={n_pos}/{len(eye_labels)}")
        except Exception as e:
            print(f"  {color}: error - {e}")

    # Also try per-SNP validation (especially rs12913832)
    print(f"\n  Per-SNP AUC (rs12913832 alone for blue vs brown):")
    if "rs12913832" in found_snps:
        row = found_snps["rs12913832"]["row"]
        snp_dose = dose[row, eye_idx_list].astype(np.float32)

        for color in ["blue", "brown"]:
            binary = np.array([1 if e == color else 0 for e in eye_labels])
            n_pos = int(np.sum(binary))
            if n_pos < 10:
                continue
            try:
                if color == "blue":
                    auc = roc_auc_score(binary, snp_dose)
                else:
                    auc = roc_auc_score(binary, -snp_dose)
                auc_best = max(auc, 1.0 - auc)
                print(f"    rs12913832 -> {color}: AUC={auc_best:.4f} (n={n_pos})")
                results.append({
                    "method": "rs12913832_only",
                    "target": color,
                    "auc": round(float(auc_best), 4),
                    "n": len(eye_labels),
                    "n_cases": n_pos,
                })
            except Exception:
                pass

    return results


def main():
    print("=" * 60)
    print("OpenSNP Validation v4 (GPU, allele-aligned)")
    print("=" * 60, flush=True)
    t_total = time.time()

    # 1. Load BIM/FAM
    print("\n[1/7] Loading BIM/FAM...", flush=True)
    bim_rsid = {}
    bim_alleles = []
    n_snps = 0
    with open(f"{HM3_PREFIX}.bim") as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            chrom, rsid, cm, pos, a1, a2 = parts[:6]
            bim_alleles.append((a1.upper(), a2.upper()))
            if rsid.startswith("rs"):
                bim_rsid[rsid] = i
            n_snps += 1

    samples = []
    with open(f"{HM3_PREFIX}.fam") as f:
        for line in f:
            samples.append(line.strip().split()[1])
    n_samples = len(samples)
    print(f"  {n_snps:,} SNPs, {n_samples} samples", flush=True)

    # 2. Read BED (now counting A1 alleles)
    print("\n[2/7] Reading BED (counting A1 alleles)...", flush=True)
    dose = read_bed(f"{HM3_PREFIX}.bed", n_samples, n_snps)

    # 3. Load weight matrix (allele-aligned v3)
    print("\n[3/7] Loading weight matrix (allele-aligned)...", flush=True)
    t_load = time.time()
    if not os.path.exists(WEIGHT_NPZ):
        # Fallback to v2 if v3 not built yet
        alt_npz = f"{BASE_DIR}/opensnp-weight-matrix.npz"
        if os.path.exists(alt_npz):
            print(f"  WARNING: v3 weight matrix not found, using v2 (NOT allele-aligned)!", flush=True)
            WEIGHT_NPZ_ACTUAL = alt_npz
        else:
            print("  ERROR: No weight matrix found!")
            sys.exit(1)
    else:
        WEIGHT_NPZ_ACTUAL = WEIGHT_NPZ

    W = sparse.load_npz(WEIGHT_NPZ_ACTUAL)
    print(f"  Shape: {W.shape}, nnz: {W.nnz:,}, loaded in {time.time()-t_load:.1f}s", flush=True)
    print(f"  Source: {WEIGHT_NPZ_ACTUAL}", flush=True)

    # Model metadata
    conn = sqlite3.connect(MEGA_DB)
    models = conn.execute("SELECT pgs_id, variant_count FROM pgs_metadata ORDER BY pgs_id").fetchall()
    trait_map = {}
    for pid, tname in conn.execute("SELECT pgs_id, trait_name FROM pgs_metadata"):
        trait_map[pid] = tname or ""
    conn.close()
    n_models = len(models)

    col_nnz = np.diff(W.indptr)
    model_coverage = {}
    for ci, (pid, vc) in enumerate(models):
        model_coverage[pid] = {
            "matched": int(col_nnz[ci]),
            "total": vc,
            "cov": round(col_nnz[ci] / max(vc, 1), 4),
        }

    # 4. GPU Scoring
    print("\n[4/7] GPU Scoring...", flush=True)
    t_score = time.time()
    scores = gpu_score(dose, W, n_models, n_samples)
    score_time = time.time() - t_score
    print(f"  {n_samples} x {n_models} scored in {score_time:.1f}s ({score_time/60:.1f}m)", flush=True)
    del W

    # 5. Ancestry
    print("\n[5/7] Assigning ancestry...", flush=True)
    ancestry = load_ancestry(PHENO_QC, samples)
    anc_counts = Counter(ancestry.values())
    print(f"  EUR: {anc_counts.get('EUR', 0)}, non-EUR: {anc_counts.get('non-EUR', 0)}, "
          f"admixed: {anc_counts.get('admixed', 0)}, unknown: {n_samples - len(ancestry)}", flush=True)

    eur_idx = [i for i, s in enumerate(samples) if ancestry.get(s) == "EUR"]
    noneur_idx = [i for i, s in enumerate(samples) if ancestry.get(s) == "non-EUR"]
    all_idx = list(range(n_samples))
    print(f"  EUR samples: {len(eur_idx)}, non-EUR: {len(noneur_idx)}", flush=True)

    # 6. Phenotypes
    print("\n[6/7] Loading phenotypes...", flush=True)
    phenotypes = {}
    with open(PHENO_FULL) as f:
        reader = csv.reader(f, delimiter=";")
        header = next(reader)
        uid_idx = header.index("user_id")
        h_idx = header.index("Height")
        hair_idx = header.index("Hair Color")
        eye_idx = header.index("Eye color")
        for row in reader:
            if len(row) <= max(h_idx, hair_idx, eye_idx):
                continue
            uid = row[uid_idx].strip()
            if uid not in phenotypes:
                phenotypes[uid] = {}
            h = parse_height(row[h_idx])
            hr = clean_hair(row[hair_idx])
            ey = clean_eye(row[eye_idx])
            if h and "height_cm" not in phenotypes[uid]:
                phenotypes[uid]["height_cm"] = h
            if hr and "hair" not in phenotypes[uid]:
                phenotypes[uid]["hair"] = hr
            if ey and "eye" not in phenotypes[uid]:
                phenotypes[uid]["eye"] = ey

    n_h = sum(1 for p in phenotypes.values() if "height_cm" in p)
    n_hr = sum(1 for p in phenotypes.values() if "hair" in p)
    n_ey = sum(1 for p in phenotypes.values() if "eye" in p)
    print(f"  Phenotypes: {len(phenotypes)} users, height={n_h}, hair={n_hr}, eye={n_ey}", flush=True)

    heights_list = [p["height_cm"] for p in phenotypes.values() if "height_cm" in p]
    if heights_list:
        print(f"  Height: mean={np.mean(heights_list):.1f}cm, std={np.std(heights_list):.1f}, "
              f"range=[{min(heights_list):.0f}, {max(heights_list):.0f}]", flush=True)
    hair_dist = Counter(p.get("hair") for p in phenotypes.values() if "hair" in p)
    print(f"  Hair: {dict(hair_dist)}", flush=True)
    eye_dist = Counter(p.get("eye") for p in phenotypes.values() if "eye" in p)
    print(f"  Eye: {dict(eye_dist)}", flush=True)

    # 7. Validation
    print("\n[7/7] Validating...", flush=True)

    # === HEIGHT ===
    print("\n" + "=" * 60)
    print("HEIGHT VALIDATION")
    print("=" * 60, flush=True)

    def get_height_arrays(idx_list):
        idxs, vals = [], []
        for i in idx_list:
            sid = samples[i]
            if sid in phenotypes and "height_cm" in phenotypes[sid]:
                idxs.append(i)
                vals.append(phenotypes[sid]["height_cm"])
        return idxs, np.array(vals) if vals else np.array([])

    h_all_idx, h_all_vals = get_height_arrays(all_idx)
    h_eur_idx, h_eur_vals = get_height_arrays(eur_idx)
    h_noneur_idx, h_noneur_vals = get_height_arrays(noneur_idx)
    print(f"Height samples: ALL={len(h_all_idx)}, EUR={len(h_eur_idx)}, non-EUR={len(h_noneur_idx)}", flush=True)

    height_kw = ["height", "stature", "body height"]
    height_results = []

    for ci, (pid, vc) in enumerate(models):
        trait = trait_map.get(pid, "")
        if not any(k in trait.lower() for k in height_kw):
            continue
        if col_nnz[ci] == 0:
            continue
        cov = model_coverage[pid]["cov"]
        result = {
            "pgs_id": pid,
            "trait": trait,
            "matched": int(col_nnz[ci]),
            "total_variants": vc,
            "coverage": cov,
        }

        if len(h_all_idx) > 30:
            sc = scores[h_all_idx, ci]
            if np.std(sc) > 1e-12:
                r, p = stats.pearsonr(sc, h_all_vals)
                result["r_all"] = round(float(r), 4)
                result["r2_all"] = round(float(r**2), 4)
                result["p_all"] = float(p)
                result["n_all"] = len(h_all_idx)

        if len(h_eur_idx) > 30:
            sc = scores[h_eur_idx, ci]
            if np.std(sc) > 1e-12:
                r, p = stats.pearsonr(sc, h_eur_vals)
                result["r_eur"] = round(float(r), 4)
                result["r2_eur"] = round(float(r**2), 4)
                result["p_eur"] = float(p)
                result["n_eur"] = len(h_eur_idx)

        if len(h_noneur_idx) > 30:
            sc = scores[h_noneur_idx, ci]
            if np.std(sc) > 1e-12:
                r, p = stats.pearsonr(sc, h_noneur_vals)
                result["r_noneur"] = round(float(r), 4)
                result["p_noneur"] = float(p)
                result["n_noneur"] = len(h_noneur_idx)

        if "r_all" in result or "r_eur" in result:
            height_results.append(result)

    height_results.sort(key=lambda x: -abs(x.get("r_eur", x.get("r_all", 0))))

    print(f"\nHeight models tested: {len(height_results)}")
    print(f"\n{'PGS ID':>14s} {'r(EUR)':>8s} {'r(ALL)':>8s} {'r(nonE)':>8s} {'cov':>6s} {'matched':>8s} {'p(EUR)':>10s}")
    print("-" * 70)
    for r in height_results[:20]:
        re_ = f"{r['r_eur']:>8.4f}" if "r_eur" in r else "     N/A"
        ra = f"{r['r_all']:>8.4f}" if "r_all" in r else "     N/A"
        rn = f"{r['r_noneur']:>8.4f}" if "r_noneur" in r else "     N/A"
        pe = f"{r['p_eur']:>10.2e}" if "p_eur" in r else "       N/A"
        print(f"{r['pgs_id']:>14s} {re_} {ra} {rn} {r['coverage']:>6.1%} {r['matched']:>8,} {pe}")

    if height_results:
        best = height_results[0]
        print(f"\nBEST HEIGHT MODEL: {best['pgs_id']}")
        print(f"  EUR:  r={best.get('r_eur', 'N/A')}, r2={best.get('r2_eur', 'N/A')}, n={best.get('n_eur', 'N/A')}")
        print(f"  ALL:  r={best.get('r_all', 'N/A')}, r2={best.get('r2_all', 'N/A')}, n={best.get('n_all', 'N/A')}")
        print(f"  Coverage: {best['coverage']:.1%} ({best['matched']:,}/{best['total_variants']:,})")

    # === HAIR COLOR ===
    print("\n" + "=" * 60)
    print("HAIR COLOR VALIDATION")
    print("=" * 60, flush=True)

    hair_kw = ["hair color", "hair colour"]
    hair_results = []
    hair_idx_list, hair_labels = [], []
    for i, sid in enumerate(samples):
        if sid in phenotypes and "hair" in phenotypes[sid]:
            hair_idx_list.append(i)
            hair_labels.append(phenotypes[sid]["hair"])
    print(f"Hair samples: {len(hair_idx_list)}, distribution: {Counter(hair_labels)}", flush=True)

    if len(hair_idx_list) > 30:
        from sklearn.metrics import roc_auc_score
        for ci, (pid, vc) in enumerate(models):
            trait = trait_map.get(pid, "")
            if not any(k in trait.lower() for k in hair_kw):
                continue
            if col_nnz[ci] == 0:
                continue
            sc = scores[hair_idx_list, ci]
            if np.std(sc) < 1e-12:
                continue
            cov = model_coverage[pid]["cov"]
            for color in ["blonde", "brown", "black", "red"]:
                binary = np.array([1 if h == color else 0 for h in hair_labels])
                n_pos = int(np.sum(binary))
                n_neg = len(binary) - n_pos
                if n_pos < 10 or n_neg < 10:
                    continue
                try:
                    auc = roc_auc_score(binary, sc)
                    if auc < 0.5:
                        auc = 1.0 - auc
                    hair_results.append({
                        "pgs_id": pid, "trait": trait, "target": color,
                        "auc": round(float(auc), 4), "n": len(hair_labels),
                        "n_cases": n_pos, "coverage": cov,
                        "matched": int(col_nnz[ci]), "total_variants": vc,
                    })
                except Exception:
                    pass

    hair_results.sort(key=lambda x: -x["auc"])
    if hair_results:
        print(f"\nHair AUC (top 20):")
        for r in hair_results[:20]:
            print(f"  {r['pgs_id']:>14s} -> {r['target']:8s} AUC={r['auc']:.4f} "
                  f"cases={r['n_cases']:>4d} cov={r['coverage']:.0%}")

    # === EYE COLOR (IrisPlex SNPs) ===
    eye_results = irisplex_eye_color(dose, bim_rsid, samples, phenotypes)

    # === BROAD TRAIT VALIDATION ===
    print("\n" + "=" * 60)
    print("BROAD TRAIT VALIDATION (all models with phenotype overlap)")
    print("=" * 60, flush=True)

    all_height_corrs = []
    if len(h_eur_idx) > 30:
        for ci, (pid, vc) in enumerate(models):
            if col_nnz[ci] == 0:
                continue
            sc = scores[h_eur_idx, ci]
            if np.std(sc) < 1e-12:
                continue
            r, p = stats.pearsonr(sc, h_eur_vals)
            if abs(r) > 0.05 and p < 0.001:
                all_height_corrs.append({
                    "pgs_id": pid,
                    "trait": trait_map.get(pid, ""),
                    "r_eur": round(float(r), 4),
                    "p_eur": float(p),
                    "matched": int(col_nnz[ci]),
                    "coverage": model_coverage[pid]["cov"],
                })
        all_height_corrs.sort(key=lambda x: -abs(x["r_eur"]))
        print(f"\nModels correlating with height (|r|>0.05, p<0.001): {len(all_height_corrs)}")
        for r in all_height_corrs[:30]:
            print(f"  {r['pgs_id']:>14s} r={r['r_eur']:>7.4f} p={r['p_eur']:.2e} "
                  f"cov={r['coverage']:.0%} | {r['trait'][:50]}")

    # === PER-ANCESTRY DISTRIBUTIONS ===
    print("\n" + "=" * 60)
    print("PER-ANCESTRY SCORE DISTRIBUTIONS")
    print("=" * 60, flush=True)

    pop_indices = {"EUR": eur_idx, "non-EUR": noneur_idx, "ALL": all_idx}
    distributions = {}
    for ci, (pid, vc) in enumerate(models):
        if col_nnz[ci] == 0:
            continue
        dist = {
            "_matched": int(col_nnz[ci]),
            "_total": vc,
            "_coverage": model_coverage[pid]["cov"],
        }
        for pop_name, pop_idx in pop_indices.items():
            if len(pop_idx) < 10:
                continue
            sc = scores[pop_idx, ci]
            if np.std(sc) < 1e-15:
                continue
            dist[pop_name] = {
                "mean": round(float(np.mean(sc)), 6),
                "std": round(float(np.std(sc)), 6),
                "min": round(float(np.min(sc)), 6),
                "max": round(float(np.max(sc)), 6),
                "p5": round(float(np.percentile(sc, 5)), 6),
                "p25": round(float(np.percentile(sc, 25)), 6),
                "p50": round(float(np.percentile(sc, 50)), 6),
                "p75": round(float(np.percentile(sc, 75)), 6),
                "p95": round(float(np.percentile(sc, 95)), 6),
                "n": len(pop_idx),
            }
        distributions[pid] = dist
    print(f"  Distributions built for {len(distributions)} models", flush=True)

    # === COVERAGE REPORT ===
    print("\n" + "=" * 60)
    print("COVERAGE REPORT")
    print("=" * 60, flush=True)

    covs = [model_coverage[pid]["cov"] for pid, vc in models if model_coverage[pid]["matched"] > 0]
    print(f"  Models with any matches: {len(covs)}/{n_models}")
    print(f"  Coverage: mean={np.mean(covs):.1%}, median={np.median(covs):.1%}")
    for t in [0.9, 0.7, 0.5, 0.3, 0.1, 0.01]:
        n = sum(1 for c in covs if c >= t)
        print(f"    >= {t:.0%}: {n} models")

    # === SAVE ===
    print("\n" + "=" * 60)
    print("SAVING RESULTS")
    print("=" * 60, flush=True)

    output = {
        "version": "v4-allele-aligned",
        "height": height_results,
        "hair_color": hair_results,
        "eye_color": eye_results,
        "broad_height_corrs": all_height_corrs,
        "distributions": distributions,
        "ancestry": {
            "EUR": len(eur_idx),
            "non-EUR": len(noneur_idx),
            "admixed": n_samples - len(eur_idx) - len(noneur_idx),
            "total": n_samples,
        },
        "coverage": {pid: model_coverage[pid] for pid, _ in models},
        "meta": {
            "n_samples": n_samples,
            "n_snps": n_snps,
            "n_models_total": n_models,
            "n_models_matched": len(covs),
            "n_height_samples_all": len(h_all_idx),
            "n_height_samples_eur": len(h_eur_idx),
            "n_hair_samples": len(hair_idx_list),
            "n_eye_samples": len(eye_results),
            "score_time_s": round(score_time, 1),
            "total_time_s": round(time.time() - t_total, 1),
            "weight_matrix": WEIGHT_NPZ_ACTUAL,
            "bed_lookup": "A1-counting [2,255,1,0]",
            "allele_aligned": "v3" in WEIGHT_NPZ_ACTUAL,
        },
    }

    with open(OUTPUT, "w") as f:
        json.dump(output, f)
    print(f"  Saved to {OUTPUT} ({os.path.getsize(OUTPUT)/1e6:.1f} MB)", flush=True)

    total_time = time.time() - t_total
    print("\n" + "=" * 60)
    print(f"DONE in {total_time:.0f}s ({total_time/60:.1f}m)")
    print("=" * 60, flush=True)


if __name__ == "__main__":
    main()
