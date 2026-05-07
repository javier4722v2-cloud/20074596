"""
Correlation analysis: Do loop genes co-express within individual SLE patients?
Uses GSE121239 (312 samples, most comprehensive dataset)
"""
import gzip, csv, json, math, os
import numpy as np
from scipy import stats as sp_stats

SM_FILE = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse121239\GSE121239_series_matrix.txt.gz'

# Known probe to gene mappings for GPL19109 (Affy U133+PM) used in GSE121239
# Based on same platform as GPL96 but extended
# Key loop gene probes from earlier full DE analysis
LOOP_PROBES = {
    # TLR pathway
    '238262_at': 'TLR7',
    '220146_at': 'TLR7',
    '222952_at': 'TLR7',
    '226906_s_at': 'TLR8',
    '220832_at': 'TLR8',
    '220832_at': 'TLR8',
    '207383_s_at': 'TLR9',
    # MYD88
    '209124_at': 'MYD88',
    # IRF7
    '208436_s_at': 'IRF7',
    # BAFF / B-cell
    '223501_at': 'TNFSF13B',  # BAFF
    '207066_at': 'TNFRSF13C',  # BAFF-R
    '207569_at': 'TNFRSF13B',  # TACI
    # IFNs
    '200887_s_at': 'STAT1',
    '209969_s_at': 'STAT1',
    '208992_s_at': 'STAT3',
    '210354_at': 'IFNG',
    # CXCL10
    '204533_at': 'CXCL10',
    # IFN signature
    '219519_s_at': 'SIGLEC1',
    '204439_at': 'IFI44L',
    '203153_at': 'IFIT1',
    '202411_at': 'IFI27',
    # TNF
    '207113_s_at': 'TNF',
    '206513_at': 'TNFSF13',   # APRIL
    # NFKB
    '209239_at': 'NFKB1',
    '211524_at': 'NFKB2',
    # Neturophils
    '206871_at': 'ELANE',
    '203949_at': 'MPO',
    # Hormones
    '203475_at': 'CYP19A1',
    '205225_at': 'ESR1',
    '211627_x_at': 'ESR1',
    # Additional
    '205841_at': 'JAK2',
    '201648_at': 'JAK1',
    '202869_at': 'OAS1',
}

# Reverse map: gene -> list of probes
GENE_PROBES = {}
for probe, gene in LOOP_PROBES.items():
    if gene not in GENE_PROBES:
        GENE_PROBES[gene] = []
    GENE_PROBES[gene].append(probe)

# Print parsing info
print("GSE121239 - Loop Gene Correlation Analysis")
print(f"Target genes: {list(GENE_PROBES.keys())}")
print(f"Total probes: {len(LOOP_PROBES)}")

# Parse series matrix
print("\nParsing series matrix (204MB)...")
sample_cols = []
gene_values = {}  # gene -> list of values per sample
probe_values = {}  # probe -> list of values per sample
header_read = False
data_started = False
n_samples = 0
probe_count = 0

with gzip.open(SM_FILE, 'rt', encoding='latin-1') as f:
    for line in f:
        line = line.rstrip('\n\r')
        
        if line.startswith('!series_matrix_table_begin'):
            data_started = True
            continue
        
        if data_started:
            if line.startswith('!series_matrix_table_end'):
                break
            
            parts = line.split('\t')
            
            if not header_read:
                # First row: sample IDs
                sample_cols = parts[1:]  # Skip ID_REF
                n_samples = len(sample_cols)
                header_read = True
                print(f"Found {n_samples} samples")
                continue
            
            probe_id = parts[0].strip()
            probe_count += 1
            
            if probe_id in LOOP_PROBES:
                gene = LOOP_PROBES[probe_id]
                values = []
                for v in parts[1:]:
                    try:
                        # Affymetrix data is typically log2 transformed
                        values.append(float(v.strip()))
                    except:
                        values.append(np.nan)
                
                if probe_id not in probe_values:
                    probe_values[probe_id] = values
            
            if probe_count % 5000 == 0:
                print(f"  Scanned {probe_count} probes...", end='\r')

print(f"\nScanned {probe_count} probes total")
print(f"Found {len(probe_values)}/{len(LOOP_PROBES)} loop probes")

# Normalize: within each sample, convert to z-score to remove batch effects
print("\nNormalizing...")
all_probe_names = list(probe_values.keys())
n_probes = len(all_probe_names)
n_samps = len(sample_cols)

# Build matrix: samples x probes
matrix = np.zeros((n_probes, n_samps))
for i, probe in enumerate(all_probe_names):
    vals = probe_values[probe]
    for j in range(min(len(vals), n_samps)):
        matrix[i, j] = vals[j]

# Replace NaN with column means
col_means = np.nanmean(matrix, axis=0)
for i in range(n_probes):
    for j in range(n_samps):
        if np.isnan(matrix[i, j]):
            matrix[i, j] = col_means[j]

# Z-score normalize each probe
for i in range(n_probes):
    mu = np.mean(matrix[i, :])
    std = np.std(matrix[i, :])
    if std > 0:
        matrix[i, :] = (matrix[i, :] - mu) / std

# Average probes per gene
gene_matrix = {}
for gene, probes in GENE_PROBES.items():
    probe_indices = [i for i, p in enumerate(all_probe_names) if p in probes]
    if probe_indices:
        avg_vals = np.mean(matrix[probe_indices, :], axis=0)
        gene_matrix[gene] = avg_vals

print(f"Genes with data: {len(gene_matrix)}")
genes_with_data = sorted(gene_matrix.keys())
print(f"Genes: {genes_with_data}")

# Compute pairwise correlations
print(f"\n=== PAIRWISE CORRELATIONS (Spearman) ===")
print(f"{'':12s}", end='')
for g2 in genes_with_data:
    print(f'{g2:>10s}', end='')
print()

corr_matrix = np.zeros((len(genes_with_data), len(genes_with_data)))
pval_matrix = np.zeros((len(genes_with_data), len(genes_with_data)))

for i, g1 in enumerate(genes_with_data):
    print(f'{g1:12s}', end='')
    for j, g2 in enumerate(genes_with_data):
        if j < i:
            print(f'      {corr_matrix[i,j]:+5.3f}', end='')
            continue
        vals1 = gene_matrix[g1]
        vals2 = gene_matrix[g2]
        r, p = sp_stats.spearmanr(vals1, vals2)
        corr_matrix[i, j] = r
        pval_matrix[i, j] = p
        if j == i:
            print(f'   1.000   ', end='')
        else:
            sig = '*' if p < 0.001 else ('·' if p < 0.05 else ' ')
            print(f'{sig}{r:+5.3f}   ', end='')
    print()

# Print significant correlations
print(f"\n=== SIGNIFICANT CORRELATIONS (p < 0.001) ===")
cutoff = 0.3
for i, g1 in enumerate(genes_with_data):
    for j, g2 in enumerate(genes_with_data):
        if j <= i:
            continue
        r = corr_matrix[i, j]
        p = pval_matrix[i, j]
        if p < 0.001:
            direction = "+++" if r > 0.5 else ("++" if r > 0.3 else "+")
            direction = direction if r > 0 else ("---" if r < -0.5 else "--" if r < -0.3 else "-")
            print(f"  {g1:10s} <-> {g2:10s}: r={r:+.3f} (p={p:.2e}) [{direction}]")

# Core LOOP CORRELATIONS (the key test)
print(f"\n=== CORE LOOP CORRELATION TEST ===")
core_loop = ['TLR7', 'MYD88', 'IRF7', 'TNFSF13B', 'CXCL10', 'STAT1', 'IFI44L']
core_idx = [i for i, g in enumerate(genes_with_data) if g in core_loop]
if len(core_idx) > 1:
    core_corr = corr_matrix[np.ix_(core_idx, core_idx)]
    flat_corr = core_corr[np.triu_indices_from(core_corr, k=1)]
    print(f"Average pairwise correlation (core loop genes): {np.mean(flat_corr):+.3f} +/- {np.std(flat_corr):.3f}")
    print(f"Median: {np.median(flat_corr):+.3f}")
    print(f"Range: [{np.min(flat_corr):+.3f}, {np.max(flat_corr):+.3f}]")
    print(f"Positive correlations: {sum(1 for c in flat_corr if c > 0)}/{len(flat_corr)}")

# IFN gene cluster correlation
print(f"\n=== IFN SIGNATURE CLUSTER ===")
ifn_genes = ['IFI44L', 'IFI27', 'IFIT1', 'SIGLEC1', 'STAT1', 'OAS1']
ifn_idx = [i for i, g in enumerate(genes_with_data) if g in ifn_genes]
if len(ifn_idx) > 1:
    ifn_corr = corr_matrix[np.ix_(ifn_idx, ifn_idx)]
    ifn_flat = ifn_corr[np.triu_indices_from(ifn_corr, k=1)]
    print(f"Average IFN cluster correlation: {np.mean(ifn_flat):+.3f} +/- {np.std(ifn_flat):.3f}")
    print(f"Positive correlations: {sum(1 for c in ifn_flat if c > 0)}/{len(ifn_flat)}")

# Check if there's a TLR7-TLR9-MYD88-BAFF-CXCL10-SIGLEC1 trajectory
print(f"\n=== PATHWAY TRAJECTORY (expected positive correlations along loop) ===")
pathway = ['TLR7', 'MYD88', 'IRF7', 'TNFSF13B', 'TNFRSF13C', 'STAT1', 'IFI44L', 'CXCL10', 'SIGLEC1']
for i, g1 in enumerate(pathway):
    if g1 not in gene_matrix:
        continue
    for j, g2 in enumerate(pathway):
        if j <= i or g2 not in gene_matrix:
            continue
        i1 = genes_with_data.index(g1)
        i2 = genes_with_data.index(g2)
        r = corr_matrix[i1, i2]
        p = pval_matrix[i1, i2]
        print(f"  {g1:12s} -> {g2:12s}: r={r:+.3f} (p={p:.2e})")

# Try to compute signature scores and correlate
print(f"\n=== SIGNATURE SCORE CORRELATIONS ===")
# IFN score: average of IFI44L, IFIT1, SIGLEC1, ISG15, OAS1
ifn_score_genes = [g for g in ['IFI44L', 'IFIT1', 'SIGLEC1', 'OAS1'] if g in gene_matrix]
if ifn_score_genes:
    ifn_score = np.mean([gene_matrix[g] for g in ifn_score_genes], axis=0)
    
# TLR score: average of TLR7 probes
tlr_score_genes = [g for g in ['TLR7'] if g in gene_matrix]
if tlr_score_genes:
    tlr_score = gene_matrix['TLR7']
    
# BAFF score
baff_genes = [g for g in ['TNFSF13B', 'TNFRSF13C'] if g in gene_matrix]
if baff_genes:
    baff_score = np.mean([gene_matrix[g] for g in baff_genes], axis=0)

if 'CXCL10' in gene_matrix:
    cxcl10 = gene_matrix['CXCL10']

print("Score correlations (Spearman):")
if tlr_score_genes and ifn_score_genes:
    r, p = sp_stats.spearmanr(tlr_score, ifn_score)
    print(f"  TLR7 score <-> IFN score: r={r:+.3f} (p={p:.2e})")
if tlr_score_genes and baff_genes:
    r, p = sp_stats.spearmanr(tlr_score, baff_score)
    print(f"  TLR7 score <-> BAFF score: r={r:+.3f} (p={p:.2e})")
if ifn_score_genes and baff_genes:
    r, p = sp_stats.spearmanr(ifn_score, baff_score)
    print(f"  IFN score  <-> BAFF score: r={r:+.3f} (p={p:.2e})")
if 'CXCL10' in gene_matrix and ifn_score_genes:
    r, p = sp_stats.spearmanr(cxcl10, ifn_score)
    print(f"  CXCL10     <-> IFN score:  r={r:+.3f} (p={p:.2e})")

# Save results
results = {
    'n_samples': n_samples,
    'genes_found': genes_with_data,
    'probes_found': len(probe_values),
    'correlation_matrix': {genes_with_data[i]: {genes_with_data[j]: float(corr_matrix[i,j]) for j in range(len(genes_with_data))} for i in range(len(genes_with_data))},
    'pval_matrix': {genes_with_data[i]: {genes_with_data[j]: float(pval_matrix[i,j]) for j in range(len(genes_with_data))} for i in range(len(genes_with_data))},
}

out = r'C:\Users\javie\.openclaw\workspace\lumus\loop_correlations.json'
with open(out, 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {out}")
