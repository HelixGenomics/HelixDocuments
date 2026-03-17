#!/usr/bin/env python3
"""Check allele alignment between PGS effect alleles and BIM A1/A2."""
import sqlite3, time

t0 = time.time()
BASE = "/workspace/pgp-pipeline"
MEGA_DB = f"{BASE}/pgs-mega.db"
BIM = f"{BASE}/opensnp-frei2024/imputed/opensnp_hm3.bim"

# Load BIM
bim_rsid = {}
bim_pos = {}
bim_alleles = []  # (a1, a2) indexed by row
with open(BIM) as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        chrom, rsid, cm, pos, a1, a2 = parts[:6]
        a1u, a2u = a1.upper(), a2.upper()
        bim_alleles.append((a1u, a2u))
        if rsid.startswith("rs"):
            bim_rsid[rsid] = i
        bim_pos[(chrom, pos)] = i

print(f"BIM loaded: {len(bim_alleles)} SNPs in {time.time()-t0:.1f}s")

# Check alignment for height models
conn = sqlite3.connect(MEGA_DB)

# Find height models
height_pgs = []
for pid, trait in conn.execute("SELECT pgs_id, trait_name FROM pgs_metadata"):
    if trait and "height" in trait.lower():
        height_pgs.append((pid, trait))

print(f"\nFound {len(height_pgs)} height models")

# Also check overall statistics across ALL models (sample 100)
all_pgs = [r[0] for r in conn.execute("SELECT pgs_id FROM pgs_metadata LIMIT 100")]

total_a1 = 0
total_a2 = 0
total_neither = 0

for pgs_id in all_pgs:
    for rsid, chrn, chrp, ea, oa, w in conn.execute(
        "SELECT rsid, chr_name, chr_position, effect_allele, other_allele, effect_weight "
        "FROM pgs_variants WHERE pgs_id = ?", (pgs_id,)):

        row_idx = None
        if rsid and rsid.startswith("rs"):
            row_idx = bim_rsid.get(rsid)
        if row_idx is None and chrn and chrp:
            row_idx = bim_pos.get((str(chrn), str(chrp)))
        if row_idx is None:
            continue

        a1, a2 = bim_alleles[row_idx]
        ea_upper = ea.upper() if ea else ""

        if ea_upper == a1:
            total_a1 += 1
        elif ea_upper == a2:
            total_a2 += 1
        else:
            total_neither += 1

total = total_a1 + total_a2 + total_neither
print(f"\n=== OVERALL (100 models) ===")
print(f"  Total matched variants: {total}")
print(f"  effect_allele == A1: {total_a1} ({total_a1/max(total,1)*100:.1f}%)")
print(f"  effect_allele == A2: {total_a2} ({total_a2/max(total,1)*100:.1f}%)")
print(f"  neither (strand/indel): {total_neither} ({total_neither/max(total,1)*100:.1f}%)")

# Check specific height models
for pid, trait in height_pgs[:3]:
    m_a1 = 0
    m_a2 = 0
    m_other = 0

    for rsid, chrn, chrp, ea, oa, w in conn.execute(
        "SELECT rsid, chr_name, chr_position, effect_allele, other_allele, effect_weight "
        "FROM pgs_variants WHERE pgs_id = ?", (pid,)):

        row_idx = None
        if rsid and rsid.startswith("rs"):
            row_idx = bim_rsid.get(rsid)
        if row_idx is None and chrn and chrp:
            row_idx = bim_pos.get((str(chrn), str(chrp)))
        if row_idx is None:
            continue

        a1, a2 = bim_alleles[row_idx]
        ea_upper = ea.upper() if ea else ""

        if ea_upper == a1:
            m_a1 += 1
        elif ea_upper == a2:
            m_a2 += 1
        else:
            m_other += 1

    m_total = m_a1 + m_a2 + m_other
    print(f"\n{pid} ({trait[:40]}):")
    print(f"  Total matched: {m_total}")
    print(f"  effect_allele == A1: {m_a1} ({m_a1/max(m_total,1)*100:.1f}%)")
    print(f"  effect_allele == A2: {m_a2} ({m_a2/max(m_total,1)*100:.1f}%)")
    print(f"  neither: {m_other} ({m_other/max(m_total,1)*100:.1f}%)")

conn.close()
print(f"\nDone in {time.time()-t0:.1f}s")
