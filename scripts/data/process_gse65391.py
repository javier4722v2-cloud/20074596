#!/usr/bin/env python3
"""
Process GSE65391 to extract expression data for target genes
and calculate fold changes between SLE and Healthy samples.
"""

import gzip
import math
import os
from collections import defaultdict

BASE = r'C:\Users\javie\.openclaw\workspace\lumus\data'

try:
    from scipy import stats as scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    # Simple t-test implementation
    def ttest_ind(a, b):
        import math as m
        n1, n2 = len(a), len(b)
        m1, m2 = sum(a)/n1, sum(b)/n2
        v1 = sum((x-m1)**2 for x in a)/(n1-1)
        v2 = sum((x-m2)**2 for x in b)/(n2-1)
        se = m.sqrt(v1/n1 + v2/n2)
        t = (m1-m2)/se if se > 0 else 0
        # Approximate p-value using normal distribution (rough)
        from scipy.stats import norm
        p = 2 * norm.sf(abs(t))
        return t, p

print("=== Step 1: Parse sample metadata ===\n")
with open(os.path.join(BASE, 'GSE65391_series_matrix.txt'), 'r', encoding='latin-1') as f:
    lines = f.readlines()

sample_gsms = []
sample_disease = []
sample_sets = []
sample_descriptions = []  # array address mapping

for line in lines:
    line = line.rstrip('\n\r')
    if line.startswith('!Sample_geo_accession'):
        parts = line.split('\t')
        sample_gsms = [p.strip().strip('"') for p in parts[1:] if p.strip()]
        print(f"Samples: {len(sample_gsms)} GSM IDs")
    elif line.startswith('!Sample_description'):
        parts = line.split('\t')
        sample_descriptions = [p.strip().strip('"') for p in parts[1:] if p.strip()]
    elif line.startswith('!Sample_characteristics_ch1') and 'disease state' in line:
        parts = line.split('\t')
        sample_disease = [p.strip().strip('"') for p in parts[1:] if p.strip()]
        h = sum(1 for s in sample_disease if 'Healthy' in s)
        s = sum(1 for s in sample_disease if 'SLE' in s)
        print(f"Disease: {h} Healthy, {s} SLE")
    elif line.startswith('!Sample_characteristics_ch1') and 'set' in line:
        parts = line.split('\t')
        sample_sets = [p.strip().strip('"') for p in parts[1:] if p.strip()]

# Extract array addresses from descriptions (format: "Sample name: 5552251039_D")
sample_to_array = {}
for i, desc in enumerate(sample_descriptions):
    if i >= len(sample_gsms):
        break
    if 'Sample name:' in desc:
        array = desc.split('Sample name:')[1].strip()
        sample_to_array[sample_gsms[i]] = array

print(f"Array mapping: {len(sample_to_array)} samples mapped")
print(f"First 5: {dict(list(sample_to_array.items())[:5])}")

# Filter out technical replicates
healthy_gsms = []
sle_gsms = []
for i, gsm in enumerate(sample_gsms):
    if i >= len(sample_disease):
        break
    is_tech = (i < len(sample_sets) and 'Technical_Replicate' in sample_sets[i])
    if not is_tech:
        if 'Healthy' in sample_disease[i]:
            healthy_gsms.append(gsm)
        elif 'SLE' in sample_disease[i]:
            sle_gsms.append(gsm)

print(f"\nAfter removing tech replicates:")
print(f"  Healthy: {len(healthy_gsms)} samples")
print(f"  SLE: {len(sle_gsms)} samples")

# Build reverse mapping: array address -> disease state
array_to_disease = {}
for i, gsm in enumerate(sample_gsms):
    if gsm in sample_to_array and i < len(sample_disease):
        array = sample_to_array[gsm]
        is_tech = (i < len(sample_sets) and 'Technical_Replicate' in sample_sets[i])
        if not is_tech:
            array_to_disease[array] = 'Healthy' if 'Healthy' in sample_disease[i] else 'SLE'

print(f"Array-disease mapping: {len(array_to_disease)} arrays")
healthy_count = sum(1 for v in array_to_disease.values() if v == 'Healthy')
sle_count = sum(1 for v in array_to_disease.values() if v == 'SLE')
print(f"  Healthy arrays: {healthy_count}, SLE arrays: {sle_count}")

# ---- Step 2: Get probe IDs for target genes ----
print("\n=== Step 2: Probe IDs for target genes ===")
target_genes = {'ESR1', 'TNF', 'TLR7', 'TNFSF13B', 'MYD88', 'CYP19A1', 'ACTB', 'GAPDH'}
probe_to_gene = {}

with open(os.path.join(BASE, 'GPL10558_full.txt'), 'r', encoding='latin-1') as f:
    in_table = False
    header = None
    for f_line in f:
        f_line = f_line.rstrip('\n\r')
        if f_line.startswith('!platform_table_begin'):
            in_table = True
            continue
        if f_line.startswith('!platform_table_end'):
            break
        if in_table:
            if header is None:
                header = f_line.split('\t')
                id_idx = header.index('ID')
                sym_idx = header.index('Symbol')
                continue
            cols = f_line.split('\t')
            if len(cols) > max(id_idx, sym_idx):
                pid = cols[id_idx]
                sym = cols[sym_idx]
                if sym in target_genes:
                    probe_to_gene[pid] = sym

target_probes = list(probe_to_gene.keys())
print(f"Target probes: {len(target_probes)}")
for p, g in probe_to_gene.items():
    print(f"  {p} -> {g}")

# ---- Step 3: Read expression data ----
print("\n=== Step 3: Reading R1 expression data ===")
r1_path = os.path.join(BASE, 'GSE65391_non-normalized_data_Illumina_HT12_V4_R1.txt.gz')

# Read expression data for target probes  
# File structure: Col 0=ID_REF, odd cols=expression values, even cols=Detection Pval
# Column headers: Col odd = array address, Col even = "Detection Pval"
probe_expr = {p: {} for p in target_probes}  # {probe: {array_addr: value}}

with gzip.open(r1_path, 'rt') as f:
    header = f.readline().strip().split('\t')
    
    # Find array address columns (odd indices starting from 1)
    array_cols = []  # (col_idx, array_address)
    for i in range(1, len(header), 2):  # Odd columns = expression data
        if i < len(header):
            arr = header[i]
            if arr in array_to_disease:  # Only if we have a mapping
                array_cols.append((i, arr))
    
    print(f"Matched array columns: {len(array_cols)} out of {len(array_to_disease)} mapped samples")
    
    lines_processed = 0
    probes_found = 0
    
    for f_line in f:
        parts = f_line.strip().split('\t')
        if len(parts) < 2:
            continue
        probe_id = parts[0]
        
        if probe_id in probe_to_gene:
            probes_found += 1
            for col_idx, arr in array_cols:
                if col_idx < len(parts):
                    try:
                        val = float(parts[col_idx])
                        probe_expr[probe_id][arr] = val
                    except ValueError:
                        pass
            if probes_found % 5 == 0:
                print(f"  Found {probes_found}/{len(target_probes)} target probes...")
        
        lines_processed += 1
    
    print(f"Processed {lines_processed} lines total")
    print(f"Found {probes_found} target probe rows")

# ---- Step 4: Calculate fold changes ----
print("\n=== Step 4: Calculating fold changes ===")

results = []
for probe_id in target_probes:
    gene = probe_to_gene[probe_id]
    expr_vals = probe_expr[probe_id]
    
    healthy_vals = [v for arr, v in expr_vals.items() if array_to_disease.get(arr) == 'Healthy']
    sle_vals = [v for arr, v in expr_vals.items() if array_to_disease.get(arr) == 'SLE']
    
    if len(healthy_vals) < 2 or len(sle_vals) < 2:
        print(f"  {probe_id} ({gene}): insufficient data (H={len(healthy_vals)}, SLE={len(sle_vals)})")
        continue
    
    mean_h = sum(healthy_vals) / len(healthy_vals)
    mean_sle = sum(sle_vals) / len(sle_vals)
    
    # The values include some negative numbers (like -7.58), which suggests
    # they may be log2 or normalized values of some kind
    # For Illumina BeadStudio output, these are typically log2 ratios
    # But negative values are unusual for log2 intensity. Let's check the range.
    all_vals = healthy_vals + sle_vals
    min_v, max_v = min(all_vals), max(all_vals)
    
    if max_v > 50:
        # Likely raw intensity values
        fc = mean_sle / mean_h if mean_h != 0 else float('inf')
        log2fc = math.log2(fc) if fc > 0 else 0
        print(f"  {probe_id} ({gene}): raw intensity range [{min_v:.1f}, {max_v:.1f}], FC={fc:.4f}")
    else:
        # Likely log2-scaled values
        fc = 2 ** (mean_sle - mean_h)
        log2fc = mean_sle - mean_h
        print(f"  {probe_id} ({gene}): log2 range [{min_v:.1f}, {max_v:.1f}], ")
        print(f"    Mean H={mean_h:.4f}, Mean SLE={mean_sle:.4f}, log2FC={log2fc:.4f}, FC={fc:.4f}")
    
    # T-test
    if HAS_SCIPY:
        t_stat, p_val = scipy_stats.ttest_ind(sle_vals, healthy_vals)
    else:
        t_stat, p_val = ttest_ind(sle_vals, healthy_vals)
    
    results.append({
        'probe': probe_id,
        'gene': gene,
        'mean_healthy': mean_h,
        'mean_sle': mean_sle,
        'fc': fc,
        'log2fc': log2fc,
        'p_val': p_val,
        'n_healthy': len(healthy_vals),
        'n_sle': len(sle_vals),
        'min_v': min_v,
        'max_v': max_v
    })
    
    print(f"    p-value = {p_val:.6e}")

# ---- Step 5: Print summary ----
print("\n\n" + "="*100)
print("RESULT SUMMARY: GSE65391 - SLE vs Healthy Fold Changes")
print("="*100)
print(f"{'Gene':<12} {'Probe':<20} {'Mean_H':<10} {'Mean_SLE':<10} {'FC':<12} {'log2FC':<10} {'p-value':<12} {'n_H':<5} {'n_SLE':<5}")
print("-"*100)
for r in sorted(results, key=lambda x: x['gene']):
    fc_display = f"{r['fc']:.4f}"
    if r['p_val'] < 0.001:
        p_display = f"{r['p_val']:.2e}"
    else:
        p_display = f"{r['p_val']:.4f}"
    print(f"{r['gene']:<12} {r['probe']:<20} {r['mean_healthy']:<10.4f} {r['mean_sle']:<10.4f} {fc_display:<12} {r['log2fc']:<10.4f} {p_display:<12} {r['n_healthy']:<5} {r['n_sle']:<5}")

print("\n\n=== Done ===")
