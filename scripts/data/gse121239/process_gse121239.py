#!/usr/bin/env python3
"""
Process GSE121239 (Toro-Dominguez 2019) series matrix to compute
fold-changes for target genes: ESR1, TNF, TLR7, BAFF, MYD88, CYP19A1
"""

import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from scipy import stats as scipy_stats

import numpy as np

DATA_DIR = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse121239'
MATRIX_FILE = f'{DATA_DIR}\\GSE121239_series_matrix.txt'

# Our target probes (standard Affy IDs for U133 Plus 2.0)
# For HT HG-U133+PM (GPL13158), probes have _PM_ suffix
TARGET_PROBES_STANDARD = {
    '207113_s_at': 'TNF',
    '223463_at': 'TLR7',
    '205225_at': 'ESR1',
    '223501_at': 'TNFSF13B (BAFF)',
    '210314_x_at': 'TNFSF13B (BAFF)',
    '220493_at': 'MYD88',
    '209124_at': 'MYD88',
    '206594_at': 'CYP19A1',
    '217398_x_at': 'GAPDH',
    '224795_x_at': 'ACTB',
}
def to_pm_id(probe_id):
    parts = probe_id.rsplit('_', 1)
    if len(parts) == 2:
        return f'{parts[0]}_PM_{parts[1]}'
    return probe_id

TARGET_PROBES = {to_pm_id(k): v for k, v in TARGET_PROBES_STANDARD.items()}

print(f"Target probes (PM platform): {list(TARGET_PROBES.keys())}")
print(f"Target genes: {set(TARGET_PROBES.values())}")

# ====== STEP 1: Parse sample metadata ======
print("\n=== Parsing sample metadata ===")

with open(MATRIX_FILE, 'r', encoding='utf-8') as f:
    lines = f.readlines()

sample_ids = []
disease_states = {}
patient_ids = {}
visit_dates = {}

for line in lines:
    if line.startswith('!Sample_title'):
        parts = line.strip().split('\t')
        for i, val in enumerate(parts[1:], start=0):
            sample_ids.append(val.strip('"'))
        print(f"Total samples: {len(sample_ids)}")
    elif line.startswith('!Sample_characteristics_ch1') and 'disease state:' in line:
        parts = line.strip().split('\t')
        for i, val in enumerate(parts[1:], start=0):
            ds = val.strip('"').replace('disease state: ', '')
            disease_states[i] = ds
        dss = [disease_states[i] for i in range(len(sample_ids))]
        c = Counter(dss)
        print(f"Disease states: {dict(c)}")
    elif line.startswith('!Sample_characteristics_ch1') and 'patient id:' in line:
        parts = line.strip().split('\t')
        for i, val in enumerate(parts[1:], start=0):
            pid = val.strip('"').replace('patient id: ', '').strip()
            patient_ids[i] = pid
        pids_unique = set(patient_ids.values())
        print(f"Unique patient IDs: {len(pids_unique)}")
    elif line.startswith('!Sample_characteristics_ch1') and 'visit date:' in line:
        parts = line.strip().split('\t')
        for i, val in enumerate(parts[1:], start=0):
            vd = val.strip('"').replace('visit date: ', '').strip()
            visit_dates[i] = vd

# Separate SLE from Healthy - FIX: check for 'Systemic Lupus Erythematosus'
sle_indices = []
healthy_indices = []
for i in range(len(sample_ids)):
    ds = disease_states.get(i, '')
    if 'SLE' in ds or 'Lupus' in ds.lower() or 'Systemic' in ds:
        sle_indices.append(i)
    elif ds == 'Healthy':
        healthy_indices.append(i)

print(f"SLE samples: {len(sle_indices)}, Healthy samples: {len(healthy_indices)}")

# Show a few SLE sample details
print("\n=== Sample details (first 25) ===")
for i in range(min(25, len(sample_ids))):
    pid = patient_ids.get(i, 'NA')
    vd = visit_dates.get(i, 'NA')
    print(f"  {i}: {sample_ids[i]:30s} | state={disease_states.get(i,'?'):35s} | patient={pid:6s} | date={vd}")

print("\n=== Sample details (SLE area) ===")
for i in range(20, min(40, len(sample_ids))):
    pid = patient_ids.get(i, 'NA')
    vd = visit_dates.get(i, 'NA')
    print(f"  {i}: {sample_ids[i]:30s} | state={disease_states.get(i,'?'):35s} | patient={pid:6s} | date={vd}")

# Count unique SLE patients from patient_ids
sle_patient_ids = set()
for i in sle_indices:
    pid = patient_ids.get(i, 'NA')
    if pid and pid != 'NA':
        sle_patient_ids.add(pid)
healthy_patient_ids = set()
for i in healthy_indices:
    pid = patient_ids.get(i, 'NA')
    if pid and pid != 'NA':
        healthy_patient_ids.add(pid)

print(f"\nUnique SLE patients (with ID): {len(sle_patient_ids)}")
print(f"Unique Healthy subjects (with ID): {len(healthy_patient_ids)}")

# Count visits per SLE patient
sle_visits_per_patient = defaultdict(list)
for i in sle_indices:
    pid = patient_ids.get(i, 'NA')
    vd = visit_dates.get(i, 'NA')
    sle_visits_per_patient[pid].append((i, vd))

visit_counts = Counter()
for pid, visits in sle_visits_per_patient.items():
    visit_counts[len(visits)] += 1
print(f"\nSLE patients by visit count: {dict(sorted(visit_counts.items()))}")

# ====== STEP 2: Parse expression matrix ======
print("\n=== Parsing expression matrix ===")

matrix_start = -1
for i, line in enumerate(lines):
    if line.startswith('!series_matrix_table_begin'):
        matrix_start = i + 1
        break

matrix_end = -1
for i in range(matrix_start, len(lines)):
    if lines[i].startswith('!series_matrix_table_end'):
        matrix_end = i
        break

print(f"Matrix rows: {matrix_end - matrix_start} (including header)")

# ====== STEP 3: Find target probes in matrix ======
print("\n=== Looking for target probes ===")

found_probes = {}
probe_data = {}

for row_idx in range(matrix_start + 1, matrix_end):
    line = lines[row_idx].strip()
    parts = line.split('\t')
    probe_id = parts[0].strip('"')
    
    if probe_id in TARGET_PROBES:
        gene = TARGET_PROBES[probe_id]
        print(f"  FOUND: {probe_id} -> {gene}")
        found_probes[probe_id] = gene
        values = {}
        for col in range(1, min(len(parts), len(sample_ids) + 1)):
            sample_idx = col - 1
            if sample_idx < len(sample_ids):
                try:
                    values[sample_idx] = float(parts[col])
                except ValueError:
                    pass
        probe_data[probe_id] = values

print(f"\nFound {len(found_probes)}/{len(TARGET_PROBES)} target probes")

# Also search for probes that might not have _PM_ suffix
print("\n=== Searching for probes without _PM_ suffix ===")
for row_idx in range(matrix_start + 1, matrix_end):
    line = lines[row_idx].strip()
    parts = line.split('\t')
    probe_id = parts[0].strip('"')
    # Match by the standard ID or any variant (e.g. 207113_s_at without PM)
    if probe_id in TARGET_PROBES_STANDARD:
        gene = TARGET_PROBES_STANDARD[probe_id]
        print(f"  FOUND (no _PM_): {probe_id} -> {gene}")
        found_probes[probe_id] = gene
        values = {}
        for col in range(1, min(len(parts), len(sample_ids) + 1)):
            sample_idx = col - 1
            if sample_idx < len(sample_ids):
                try:
                    values[sample_idx] = float(parts[col])
                except ValueError:
                    pass
        probe_data[probe_id] = values
    # Also try _at only (some probes might not have _s_, _x_ prefix variations)
    base_id = probe_id.replace('_PM_', '_')
    if base_id in TARGET_PROBES_STANDARD and probe_id not in found_probes and probe_id not in [p for p in found_probes]:
        if base_id != probe_id:  # only if it was a PM variant
            gene = TARGET_PROBES_STANDARD[base_id]
            found_probes[probe_id] = gene
            values = {}
            for col in range(1, min(len(parts), len(sample_ids) + 1)):
                sample_idx = col - 1
                if sample_idx < len(sample_ids):
                    try:
                        values[sample_idx] = float(parts[col])
                    except ValueError:
                        pass
            probe_data[probe_id] = values

# ====== STEP 4: Handle longitudinal data ======
print("\n=== Handling longitudinal data ===")

# For each patient, find their first visit
sle_patient_samples = defaultdict(list)
healthy_patient_samples = defaultdict(list)

for i in range(len(sample_ids)):
    ds = disease_states.get(i, '')
    pid = patient_ids.get(i, 'NA')
    vd = visit_dates.get(i, 'NA')
    if 'SLE' in ds or 'Lupus' in ds.lower() or 'Systemic' in ds:
        sle_patient_samples[pid].append((i, vd))
    elif ds == 'Healthy':
        healthy_patient_samples[pid].append((i, vd))

print(f"SLE patient groups: {len(sle_patient_samples)}")
print(f"Healthy subject groups: {len(healthy_patient_samples)}")

def parse_date(date_str):
    if date_str == 'NA' or not date_str:
        return None
    # Detect format: YYYY-MM-DD or DD/MM/YYYY
    if '/' in date_str:
        parts = date_str.split('/')
        return (int(parts[2]), int(parts[1]), int(parts[0]))
    elif '-' in date_str:
        parts = date_str.split('-')
        return (int(parts[0]), int(parts[1]), int(parts[2]))
    return None

def get_first_visit(samples):
    """Return sample index of first visit"""
    if len(samples) == 1:
        return samples[0][0]
    
    dated_samples = []
    for idx, date_str in samples:
        d = parse_date(date_str)
        if d:
            dated_samples.append((d, idx))
    
    if dated_samples:
        dated_samples.sort(key=lambda x: x[0])
        return dated_samples[0][1]
    
    # No dates available, use first in list
    return samples[0][0]

# Get first visit for each SLE patient
sle_first_visit = []
for pid, samples in sle_patient_samples.items():
    idx = get_first_visit(samples)
    sle_first_visit.append(idx)

# Healthy: if multiple visits (shouldn't be, but just in case), take first
healthy_unique = []
for pid, samples in healthy_patient_samples.items():
    if pid == 'NA' or len(samples) == 1:
        for idx, _ in samples:
            healthy_unique.append(idx)
    else:
        idx = get_first_visit(samples)
        healthy_unique.append(idx)

# Remove duplicates
sle_first_visit = sorted(set(sle_first_visit))
healthy_unique = sorted(set(healthy_unique))

print(f"\nAfter deduplication (first visit only):")
print(f"  SLE samples: {len(sle_first_visit)}")
print(f"  Healthy samples: {len(healthy_unique)}")

# Show which samples were selected
print("\nFirst-visit SLE samples:")
for i in sle_first_visit:
    pid = patient_ids.get(i, 'NA')
    vd = visit_dates.get(i, 'NA')
    print(f"  {i}: {sample_ids[i]:30s} | patient={pid:6s} | date={vd}")

print("\nHealthy samples:")
for i in healthy_unique:
    pid = patient_ids.get(i, 'NA')
    print(f"  {i}: {sample_ids[i]:30s} | patient={pid:6s}")

# ====== STEP 5: Calculate fold changes ======
print("\n=== Calculating fold changes ===")

results = []

for probe_id, gene in sorted(found_probes.items(), key=lambda x: x[1]):
    if probe_id not in probe_data:
        print(f"  WARNING: {probe_id} ({gene}) not found in expression data!")
        continue
    
    vals = probe_data[probe_id]
    
    sle_vals = [vals[i] for i in sle_first_visit if i in vals]
    healthy_vals = [vals[i] for i in healthy_unique if i in vals]
    
    if len(sle_vals) < 2 or len(healthy_vals) < 2:
        print(f"  WARNING: {gene} ({probe_id}): too few samples (SLE={len(sle_vals)}, Healthy={len(healthy_vals)})")
        continue
    
    sleep_mean = np.mean(sle_vals)
    healthy_mean = np.mean(healthy_vals)
    
    fc_log2 = sleep_mean - healthy_mean
    fc_linear = 2 ** fc_log2
    
    t_stat, p_val = scipy_stats.ttest_ind(sle_vals, healthy_vals, equal_var=False)
    
    if fc_log2 > 0:
        direction = "UP in SLE"
    else:
        direction = "DOWN in SLE"
    
    results.append({
        'gene': gene,
        'probe': probe_id,
        'sle_mean': sleep_mean,
        'healthy_mean': healthy_mean,
        'fc_log2': fc_log2,
        'fc_linear': fc_linear,
        'p_value': p_val,
        'direction': direction,
        'n_sle': len(sle_vals),
        'n_healthy': len(healthy_vals),
    })
    
    print(f"  {gene:20s} ({probe_id:25s}): SLE={sleep_mean:.3f} Ctrl={healthy_mean:.3f} "
          f"FC={fc_linear:.4f} (log2={fc_log2:+.4f}) p={p_val:.4e} "
          f"{direction} (n_SLE={len(sle_vals)}, n_Ctrl={len(healthy_vals)})")

# ====== STEP 6: Output summary ======
print("\n\n" + "="*90)
print("RESULT SUMMARY")
print("="*90)
header = f"{'Gene':22s} {'Probe':25s} {'SLE_mean':10s} {'Ctrl_mean':10s} {'FC_lin':10s} {'FC_log2':10s} {'p-val':12s} {'Dir':12s}"
print(header)
print("-"*111)
for r in sorted(results, key=lambda x: x['gene']):
    p_str = f"{r['p_value']:.4e}" if r['p_value'] > 1e-100 else "<1e-100"
    print(f"{r['gene']:22s} {r['probe']:25s} {r['sle_mean']:10.4f} {r['healthy_mean']:10.4f} "
          f"{r['fc_linear']:10.4f} {r['fc_log2']:10.4f} {p_str:12s} {r['direction']:12s}")

# Save results to JSON
import json
output = []
for r in results:
    output.append({
        'gene': r['gene'],
        'probe': r['probe'],
        'sle_mean': round(r['sle_mean'], 4),
        'healthy_mean': round(r['healthy_mean'], 4),
        'fc_linear': round(r['fc_linear'], 4),
        'fc_log2': round(r['fc_log2'], 4),
        'p_value': r['p_value'],
        'direction': r['direction'],
        'n_sle': r['n_sle'],
        'n_healthy': r['n_healthy'],
    })

with open(f'{DATA_DIR}\\gse121239_results.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {DATA_DIR}\\gse121239_results.json")
