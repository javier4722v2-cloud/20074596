"""Refined CLE analysis - correct FC direction, all comparisons"""
import pandas as pd
import numpy as np
import os
from scipy import stats
from collections import defaultdict
import json

DATA_DIR = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data'
OUT_DIR = os.path.join(DATA_DIR, 'results')
os.makedirs(OUT_DIR, exist_ok=True)

# Load probe-to-gene mapping
probe_to_gene = {}
with open(os.path.join(DATA_DIR, 'probe_to_gene.txt'), encoding='utf-8') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            probe_to_gene[parts[0]] = parts[1]

def collapse_probes_to_genes(df, mapping):
    from collections import defaultdict
    gene_probes = defaultdict(list)
    for probe in df.index:
        if probe in mapping:
            gene_probes[mapping[probe]].append(probe)
    gene_expr = {}
    for gene, probes in gene_probes.items():
        if len(probes) == 1:
            gene_expr[gene] = df.loc[probes[0]]
        else:
            gene_expr[gene] = df.loc[probes].median()
    return pd.DataFrame(gene_expr).T

def classify_sample(title):
    t = title.lower()
    if 'chronic cutaneous' in t: return 'CCLE'
    elif 'subacute cutaneous' in t: return 'SCLE'
    elif 'acute cutaneous' in t: return 'ACLE'
    elif 'control' in t or 'normal' in t: return 'Control'
    elif 'psoriasis' in t: return 'Psoriasis'
    elif 'lupus nephritis' in t: return 'LN'
    elif 'kidney control' in t: return 'KidneyControl'
    return 'Other'

def de_analysis(df, group_label, samples_case, samples_ctrl):
    """Calculate FC = case/control, log2FC positive = UP in case"""
    results = []
    for gene in df.index:
        vals_case = df.loc[gene, samples_case].values.astype(float)
        vals_ctrl = df.loc[gene, samples_ctrl].values.astype(float)
        vc = vals_case[~np.isnan(vals_case)]
        vt = vals_ctrl[~np.isnan(vals_ctrl)]
        if len(vc) < 2 or len(vt) < 2:
            continue
        if np.mean(vc) <= 0 or np.mean(vt) <= 0:
            continue
        fc = np.mean(vc) / np.mean(vt)
        l2fc = np.log2(fc) if fc > 0 else 0
        t_stat, p_val = stats.ttest_ind(vc, vt)
        results.append({'gene': gene, 'mean_case': np.mean(vc), 'mean_ctrl': np.mean(vt),
                        'fold_change': fc, 'log2FC': l2fc, 't_stat': t_stat, 'p_value': p_val})
    
    out = pd.DataFrame(results)
    if len(out) > 0:
        out['abs_t'] = out['t_stat'].abs()
        out = out.sort_values('abs_t', ascending=False)
    return out

def print_top_de(de, sle_genes, label, n=15):
    print(f"\n  === {label} ===")
    print(f"  Top {n} DE genes (case/control FC, by |t-stat|):")
    for _, r in de.head(n).iterrows():
        d = '↑' if r['fold_change'] > 1 else '↓'
        print(f"    {r['gene']:15s}  FC={r['fold_change']:.2f}  log2={r['log2FC']:+.2f}  p={r['p_value']:.1e}  {d}")

    sig = de[de['p_value'] < 0.05]
    up = sig[sig['fold_change'] > 1]
    dn = sig[sig['fold_change'] < 1]
    print(f"  Significant (p<0.05): {len(sig)} genes ({len(up)} UP, {len(dn)} DOWN)")
    
    print(f"  SLE Loop Genes:")
    for gene in sle_genes:
        rw = de[de['gene'] == gene]
        if len(rw) > 0:
            r = rw.iloc[0]
            fl = '↑' if r['fold_change'] > 1 else '↓'
            p = '*' if r['p_value'] < 0.05 else ''
            print(f"    {gene:15s}  FC={r['fold_change']:.2f}  log2={r['log2FC']:+.2f}  p={r['p_value']:.1e}  {fl}{p}")
        else:
            print(f"    {gene:15s}  NOT IN DATA")

# ======= MAIN =======
sle_genes = ['TLR7','TLR8','MYD88','TNFSF13B','IRF7','CXCL10','ELANE','TNF',
             'STAT1','STAT3','IL6','IL1B','TLR9','CYP19A1','ESR1','IRAK4','TICAM1']

# --- Dataset 1: GSE109248 ---
print("="*65)
print("DATASET 1: CLE SKIN BIOPSY (GSE109248)")
print("="*65)

df109, titles109, geo109 = None, [], []
with open(os.path.join(DATA_DIR, 'GSE109248_series_matrix.txt'), encoding='utf-8', errors='replace') as f:
    text = f.read()
lines = text.split('\n')
for line in lines:
    if line.startswith('!Sample_title'): titles109 = [s.strip().strip('"') for s in line.split('\t')[1:]]
    if line.startswith('!Sample_geo_accession'): geo109 = [s.strip().strip('"') for s in line.split('\t')[1:]]
start = next(i for i,l in enumerate(lines) if 'series_matrix_table_begin' in l) + 1
data_lines = []
for i in range(start, len(lines)):
    if 'series_matrix_table_end' in lines[i]: break
    data_lines.append(lines[i])
from io import StringIO
df109 = pd.read_csv(StringIO('\n'.join(data_lines)), sep='\t', index_col=0)
df109.columns = [c.strip().strip('"') for c in df109.columns]

genes109 = collapse_probes_to_genes(df109, probe_to_gene)
groups109 = [classify_sample(t) for t in titles109]
sample_groups = defaultdict(list)
for gid, grp in zip(geo109, groups109):
    sample_groups[grp].append(gid)

print(f"Probes: {df109.shape[0]} → Genes: {genes109.shape[0]}")
for g in ['CCLE','SCLE','ACLE','Control','Psoriasis']:
    print(f"  {g}: {len(sample_groups.get(g,[]))} samples")

# Comparison 1: CLE (all) vs Control
a = sample_groups['CCLE']+sample_groups['SCLE']+sample_groups['ACLE']
b = sample_groups['Control']
if a and b:
    de1 = de_analysis(genes109, 'CLE_vs_Control', a, b)
    de1.to_csv(os.path.join(OUT_DIR, '109_CLE_vs_Control.csv'), index=False)
    print_top_de(de1, sle_genes, 'CLE vs Control (skin)')

# Comparison 2: CLE vs Psoriasis (both skin autoimmune)
a = sample_groups['CCLE']+sample_groups['SCLE']+sample_groups['ACLE']
b = sample_groups['Psoriasis']
if a and b:
    de2 = de_analysis(genes109, 'CLE_vs_Psoriasis', a, b)
    de2.to_csv(os.path.join(OUT_DIR, '109_CLE_vs_Psoriasis.csv'), index=False)
    print_top_de(de2, sle_genes, 'CLE vs Psoriasis (both skin autoimmune)')

# Comparison 3: CCLE vs SCLE
a, b = sample_groups['CCLE'], sample_groups['SCLE']
if a and b:
    de3 = de_analysis(genes109, 'CCLE_vs_SCLE', a, b)
    de3.to_csv(os.path.join(OUT_DIR, '109_CCLE_vs_SCLE.csv'), index=False)
    print_top_de(de3, sle_genes, 'CCLE vs SCLE')

# --- Dataset 2: GSE167923 (CLE blood, no controls) ---
print("\n" + "="*65)
print("DATASET 2: CLE WHOLE BLOOD RNA-SEQ (GSE167923)")
print("="*65)

expr = pd.read_csv(os.path.join(DATA_DIR, 'GSE167923_expression.csv'))
gene_col = expr.columns[0]
expr = expr.set_index(gene_col)
cols = list(expr.columns)
print(f"Samples: {len(cols)} (no controls - CLE only)")

# From the paper: 62 CLE patients, various subtypes
# Without control samples, we can look at:
# - Range of expression across CLE patients
# - Correlation structure

# Show expression range for SLE loop genes
print(f"\n  Expression range for SLE loop genes (across {len(cols)} CLE patients):")
for gene in sle_genes:
    if gene in expr.index:
        vals = expr.loc[gene].values.astype(float)
        print(f"    {gene:15s}  mean={np.mean(vals):.1f}  range=[{np.min(vals):.1f}, {np.max(vals):.1f}]  CV={np.std(vals)/np.mean(vals):.2f}")

# Look for internal clusters - high variance genes
variances = expr.var(axis=1).sort_values(ascending=False)
print(f"\n  Top 10 most variable genes across CLE patients:")
for g in variances.head(10).index:
    print(f"    {g:15s}  var={variances[g]:.2f}")

# --- Dataset 3: GSE112943 (CLE skin + LN kidney) ---
print("\n" + "="*65)
print("DATASET 3: CLE SKIN + LUPUS NEPHRITIS (GSE112943)")
print("="*65)

df112, titles112, geo112 = None, [], []
with open(os.path.join(DATA_DIR, 'GSE112943_series_matrix.txt'), encoding='utf-8', errors='replace') as f:
    text = f.read()
lines = text.split('\n')
for line in lines:
    if line.startswith('!Sample_title'): titles112 = [s.strip().strip('"') for s in line.split('\t')[1:]]
    if line.startswith('!Sample_geo_accession'): geo112 = [s.strip().strip('"') for s in line.split('\t')[1:]]
start = next(i for i,l in enumerate(lines) if 'series_matrix_table_begin' in l) + 1
data_lines = []
for i in range(start, len(lines)):
    if 'series_matrix_table_end' in lines[i]: break
    data_lines.append(lines[i])
df112 = pd.read_csv(StringIO('\n'.join(data_lines)), sep='\t', index_col=0)
df112.columns = [c.strip().strip('"') for c in df112.columns]

genes112 = collapse_probes_to_genes(df112, probe_to_gene)
groups112 = [classify_sample(t) for t in titles112]
sample_groups112 = defaultdict(list)
for gid, grp in zip(geo112, groups112):
    sample_groups112[grp].append(gid)

print(f"Probes: {df112.shape[0]} → Genes: {genes112.shape[0]}")
for g in ['CCLE','SCLE','Control','LN','KidneyControl']:
    print(f"  {g}: {len(sample_groups112.get(g,[]))} samples")

# CLE skin vs control skin (replicate from GSE109248)
a = sample_groups112['CCLE']+sample_groups112['SCLE']
b = sample_groups112['Control']
if a and b:
    de4 = de_analysis(genes112, '112_CLE_vs_Control_skin', a, b)
    de4.to_csv(os.path.join(OUT_DIR, '112_CLE_vs_Control_skin.csv'), index=False)
    print_top_de(de4, sle_genes, 'CLE vs Control (skin, replicate)')

# LN vs Control kidney
a = sample_groups112['LN']
b = sample_groups112['KidneyControl']
if a and b:
    de5 = de_analysis(genes112, '112_LN_vs_KidneyControl', a, b)
    de5.to_csv(os.path.join(OUT_DIR, '112_LN_vs_KidneyControl.csv'), index=False)
    print_top_de(de5, sle_genes, 'LN vs Control (kidney)')

# --- Cross-dataset: Gather SLE loop gene behavior ---
print("\n" + "="*65)
print("CROSS-DATASET SUMMARY: SLE LOOP GENES")
print("="*65)
print(f"\n{'Gene':15s} {'SLE Blood':12s} {'CLE Skin':12s} {'CLE Skin2':12s} {'CLE Blood':12s} {'Psor Skin':12s} {'LN Kidney':12s}")
print("-"*85)

for gene in sle_genes:
    vals = {}
    for tag, de, direction in [
        ('SLE Blood', None, '↑'):],
    # We don't have SLE blood DE here - would need to match with our existing datasets
    # For now, show CLE comparisons
        pass
    
# Actually let me just compile from the saved results
print("\n(To add SLE blood comparison, need to integrate with our existing SLE datasets)")

print("\n\nFiles saved to:", OUT_DIR)
print("Done!")
