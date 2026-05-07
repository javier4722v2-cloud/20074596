"""
Analyze GSE81622: Extract FCs for target genes (SLE vs Control)
"""
import re
import csv
import math
from collections import defaultdict
from scipy import stats as scipy_stats

# Probe to Gene mapping (GPL10558)
probes = {
    'ILMN_1728106': 'TNF',
    'ILMN_1677827': 'TLR7',
    'ILMN_1678535': 'ESR1',
    'ILMN_1758418': 'TNFSF13B(BAFF)',
    'ILMN_2066858': 'TNFSF13B(BAFF)',
    'ILMN_1738523': 'MYD88',
    'ILMN_1699139': 'CYP19A1',
    'ILMN_1343295': 'GAPDH',
    'ILMN_1802252': 'GAPDH',
    'ILMN_2038778': 'GAPDH',
    'ILMN_1777296': 'ACTB',
    'ILMN_2152131': 'ACTB',
}

target_probes = set(probes.keys())

# Build description -> condition mapping from metadata
desc_to_condition = {}
with open(r'C:\Users\javie\.openclaw\workspace\lumus\data\gse81622\GSE81622_metadata.txt', 'r', encoding='latin-1') as f:
    data = f.read()

samples_raw = data.split('^SAMPLE = ')[1:]
for s in samples_raw:
    title = ''
    desc_lines = []
    for line in s.split('\n'):
        if line.startswith('!Sample_title = '):
            title = line.split('= ')[1].strip()
        elif line.startswith('!Sample_description = '):
            desc_lines.append(line.split('= ')[1].strip())
    
    # Determine condition
    if 'normal control' in title.lower():
        condition = 'Control'
    elif 'SLE' in title:
        condition = 'SLE'
    else:
        continue
    
    for dl in desc_lines:
        if 'replicate' not in dl and dl != 'NONE':
            desc_to_condition[dl] = condition

print(f"Sample mapping: {len(desc_to_condition)} samples")
ctrl_count = sum(1 for v in desc_to_condition.values() if v == 'Control')
sle_count = sum(1 for v in desc_to_condition.values() if v == 'SLE')
print(f"Controls: {ctrl_count}, SLE: {sle_count}")

# Read data file and extract probe values
# Structure: ID_REF, then groups of 4 columns per sample: AVG_Signal, Avg_NBEADS, BEAD_STDERR, Detection Pval
probe_data = defaultdict(lambda: defaultdict(list))  # {gene: {condition: [values]}}

with open(r'C:\Users\javie\.openclaw\workspace\lumus\data\gse81622\GSE81622_non-normalized.txt', 'r', encoding='latin-1') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    
    # Build column mapping: sample description -> (avg_signal_col, detection_pval_col)
    sample_cols = {}
    for i in range(1, len(header), 4):
        desc = header[i].replace('.AVG_Signal', '')
        sample_cols[desc] = {
            'avg_col': i,
            'pval_col': i + 3
        }
    
    # Check which sample descriptions we have
    found_mapping = sum(1 for d in sample_cols if d in desc_to_condition)
    print(f"Sample descriptions mapped to conditions: {found_mapping} / {len(sample_cols)}")
    unmapped = [d for d in sample_cols if d not in desc_to_condition]
    if unmapped:
        print(f"Unmapped samples: {unmapped}")
    
    # Process each row
    row_count = 0
    for row in reader:
        if not row:
            continue
        probe_id = row[0]
        if probe_id in target_probes:
            gene = probes[probe_id]
            row_count += 1
            
            for desc, cols in sample_cols.items():
                condition = desc_to_condition.get(desc)
                if condition:
                    try:
                        val = float(row[cols['avg_col']])
                        probe_data[probe_id][condition].append(val)
                    except (ValueError, IndexError):
                        pass

print(f"\nFound {row_count} probe rows matching targets")

# Aggregate probes per gene (average across probes for same gene)
gene_data = defaultdict(lambda: defaultdict(list))  # {gene: {condition: [values]}}

for probe_id, cond_data in probe_data.items():
    gene = probes[probe_id]
    for condition, values in cond_data.items():
        gene_data[gene][condition].extend(values)

# Calculate statistics
print("\n==============================")
print("RESULTS: GSE81622")
print("==============================\n")
print(f"{'Gene':<20} {'SLE_mean':<12} {'Ctrl_mean':<12} {'FC':<12} {'p_val':<12} {'log2FC':<12}")
print("-" * 80)

results = []
for gene in ['TNF', 'TLR7', 'ESR1', 'TNFSF13B(BAFF)', 'MYD88', 'CYP19A1', 'GAPDH', 'ACTB']:
    if gene not in gene_data:
        print(f"{gene:<20} {'N/A':<12} {'N/A':<12} {'N/A':<12} {'N/A':<12} {'N/A':<12}")
        continue
    
    sle_vals = gene_data[gene].get('SLE', [])
    ctrl_vals = gene_data[gene].get('Control', [])
    
    if not sle_vals or not ctrl_vals:
        print(f"{gene:<20} {'NO DATA':<12}")
        continue
    
    sle_mean = sum(sle_vals) / len(sle_vals)
    ctrl_mean = sum(ctrl_vals) / len(ctrl_vals)
    
    # Check if data is in log2 space (typical Avg_Signal values are 50-20000, so linear)
    fc = sle_mean / ctrl_mean if ctrl_mean > 0 else float('inf')
    log2fc = math.log2(fc) if fc > 0 else float('-inf')
    
    # Welch's t-test
    t_stat, p_val = scipy_stats.ttest_ind(sle_vals, ctrl_vals, equal_var=False)
    
    results.append({
        'gene': gene,
        'sle_mean': sle_mean,
        'ctrl_mean': ctrl_mean,
        'fc': fc,
        'log2fc': log2fc,
        'p_val': p_val,
        'n_sle': len(sle_vals),
        'n_ctrl': len(ctrl_vals)
    })
    
    print(f"{gene:<20} {sle_mean:<12.2f} {ctrl_mean:<12.2f} {fc:<12.4f} {p_val:<12.6f} {log2fc:<12.4f}")

print("\n\n--- Detailed per-probe values ---")
for probe_id in sorted(probe_data.keys()):
    gene = probes[probe_id]
    sle_vals = probe_data[probe_id].get('SLE', [])
    ctrl_vals = probe_data[probe_id].get('Control', [])
    
    if sle_vals and ctrl_vals:
        sle_mean = sum(sle_vals) / len(sle_vals)
        ctrl_mean = sum(ctrl_vals) / len(ctrl_vals)
        fc = sle_mean / ctrl_mean if ctrl_mean > 0 else float('inf')
        t, p = scipy_stats.ttest_ind(sle_vals, ctrl_vals, equal_var=False)
        print(f"{probe_id} ({gene}): SLE={sle_mean:.2f}(n={len(sle_vals)}) Ctrl={ctrl_mean:.2f}(n={len(ctrl_vals)}) FC={fc:.4f} p={p:.6f}")

print("\n\n--- Comparison with GSE50772 ---")
print(f"{'Gene':<20} {'GSE81622_FC':<12} {'GSE50772_FC':<12} {'Match':<12}")
print("-" * 56)
# GSE50772 values (from the task description)
gse50772_fc = {
    'TNF': 4.3,
    'ESR1': 0.38,
    'TNFSF13B(BAFF)': 1.67,
}

for r in results:
    gene = r['gene']
    if gene in gse50772_fc:
        ref_fc = gse50772_fc[gene]
        match = '✓' if (r['fc'] > 1 and ref_fc > 1) or (r['fc'] < 1 and ref_fc < 1) else '✗'
        print(f"{gene:<20} {r['fc']:<12.4f} {ref_fc:<12.2f} {match:<12}")
