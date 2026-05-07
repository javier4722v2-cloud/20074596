"""Correlation analysis: CLE vs SLE loop gene correlations"""
import pandas as pd
import numpy as np
import os
from scipy import stats

DATA_DIR = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data'

sle_loop = ['TLR7','TLR8','MYD88','TNFSF13B','IRF7','CXCL10','ELANE',
            'STAT1','STAT3','IL1B','TLR9','IRAK4','TICAM1','IL6','TNF']
ifn_genes = ['IFI44','ISG15','IFI44L','IFIT3','OAS1','OASL','IFIT1','IFIH1','RSAD2','HERC5']
all_genes = sle_loop + ifn_genes

def corr_matrix(df_genes, sample_cols, label):
    """Pearson correlation between loop genes across samples"""
    vals = {}
    for g in all_genes:
        if g in df_genes.index:
            v = df_genes.loc[g, sample_cols].values.astype(float)
            v = v[~np.isnan(v)]
            if len(v) >= len(sample_cols)*0.5:  # at least 50% valid
                vals[g] = v
    
    print(f"\n{'='*60}")
    print(f"  {label} ({len(sample_cols)} samples)")
    print(f"{'='*60}")
    
    genes = list(vals.keys())
    if len(genes) < 2:
        print("  Too few genes for correlation")
        return
    
    # Build correlation matrix
    n = len(genes)
    corr_mat = pd.DataFrame(np.zeros((n, n)), index=genes, columns=genes)
    pval_mat = pd.DataFrame(np.ones((n, n)), index=genes, columns=genes)
    
    for i, g1 in enumerate(genes):
        for j, g2 in enumerate(genes):
            if i >= j:
                continue
            v1, v2 = vals[g1], vals[g2]
            min_len = min(len(v1), len(v2))
            if min_len < 5:
                continue
            # Align by shortest
            r, p = stats.pearsonr(v1[:min_len], v2[:min_len])
            corr_mat.loc[g1, g2] = r
            corr_mat.loc[g2, g1] = r
            pval_mat.loc[g1, g2] = p
            pval_mat.loc[g2, g1] = p
    
    # Print key correlations (the loop: TLR7→MYD88→IRF7→CXCL10, BAFF)
    print(f"\n  TLR7 correlations:")
    if 'TLR7' in genes:
        for g in genes:
            if g != 'TLR7':
                r = corr_mat.loc['TLR7', g]
                p = pval_mat.loc['TLR7', g]
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                print(f"    TLR7 vs {g:12s}: r={r:+.3f}  p={p:.1e} {sig}")
    
    print(f"\n  IRF7 correlations (key loop hub):")
    if 'IRF7' in genes:
        for g in genes:
            if g != 'IRF7':
                r = corr_mat.loc['IRF7', g]
                p = pval_mat.loc['IRF7', g]
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                print(f"    IRF7 vs {g:12s}: r={r:+.3f}  p={p:.1e} {sig}")
    
    print(f"\n  BAFF/TNFSF13B correlations:")
    if 'TNFSF13B' in genes:
        for g in genes:
            if g != 'TNFSF13B':
                r = corr_mat.loc['TNFSF13B', g]
                p = pval_mat.loc['TNFSF13B', g]
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                print(f"    BAFF vs {g:12s}: r={r:+.3f}  p={p:.1e} {sig}")
    
    print(f"\n  CXCL10 correlations (output marker):")
    if 'CXCL10' in genes:
        for g in genes:
            if g != 'CXCL10':
                r = corr_mat.loc['CXCL10', g]
                p = pval_mat.loc['CXCL10', g]
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                print(f"    CXCL10 vs {g:12s}: r={r:+.3f}  p={p:.1e} {sig}")
    
    # Heatmap summary
    print(f"\n  Correlation heatmap (all vs all):")
    print("       ", "".join(f"{g:>7s}" for g in genes))
    for g1 in genes:
        row = f"  {g1:>7s}"
        for g2 in genes:
            r = corr_mat.loc[g1, g2]
            if g1 == g2:
                row += "   1.00"
            else:
                row += f" {r:+.2f} " if r >= 0 else f"{r:+.2f} "
        print(row)
    
    return corr_mat

# ====== 1. GSE109248: CLE skin samples ======
print("\n\n1. CLE SKIN - GSE109248")
print("="*60)

# Load and parse
df, titles, geo = None, [], []
with open(os.path.join(DATA_DIR, 'GSE109248_series_matrix.txt'), encoding='utf-8', errors='replace') as f:
    text = f.read()
lines = text.split('\n')
for l in lines:
    if l.startswith('!Sample_title'): titles = [s.strip().strip('"') for s in l.split('\t')[1:]]
    if l.startswith('!Sample_geo_accession'): geo = [s.strip().strip('"') for s in l.split('\t')[1:]]
start = next(i for i,l in enumerate(lines) if 'series_matrix_table_begin' in l) + 1
dl = []
for i in range(start, len(lines)):
    if 'series_matrix_table_end' in lines[i]: break
    dl.append(lines[i])
from io import StringIO
df = pd.read_csv(StringIO('\n'.join(dl)), sep='\t', index_col=0)
df.columns = [c.strip().strip('"') for c in df.columns]

# Load probe-to-gene mapping
probe_to_gene = {}
with open(os.path.join(DATA_DIR, 'probe_to_gene.txt'), encoding='utf-8') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            probe_to_gene[parts[0]] = parts[1]

from collections import defaultdict
gp = defaultdict(list)
for p in df.index:
    if p in probe_to_gene:
        gp[probe_to_gene[p]].append(p)
gene_expr = {}
for g, ps in gp.items():
    gene_expr[g] = df.loc[ps[0]] if len(ps) == 1 else df.loc[ps].median()
genes109 = pd.DataFrame(gene_expr).T

# Classify samples
def cls(title):
    t = title.lower()
    if 'chronic' in t: return 'CCLE'
    elif 'subacute' in t: return 'SCLE'
    elif 'acute' in t: return 'ACLE'
    elif 'psoriasis' in t: return 'Psoriasis'
    elif 'control' in t or 'normal' in t: return 'Control'
    return 'Other'

groups = [cls(t) for t in titles]

# Correlation in CLE skin samples (all subtypes)
cle_samples = [geo[i] for i, g in enumerate(groups) if g in ['CCLE','SCLE','ACLE']]
print(f"\nCLE skin samples: {len(cle_samples)}")
corr_matrix(genes109, cle_samples, "CLE skin - GSE109248")

# Correlation in control skin (baseline)
ctrl_samples = [geo[i] for i, g in enumerate(groups) if g == 'Control']
print(f"\n\nControl skin samples: {len(ctrl_samples)}")
corr_matrix(genes109, ctrl_samples, "Control skin - GSE109248")

# Correlation in psoriasis skin (other autoimmune skin)
pso_samples = [geo[i] for i, g in enumerate(groups) if g == 'Psoriasis']
print(f"\n\nPsoriasis skin samples: {len(pso_samples)}")
corr_matrix(genes109, pso_samples, "Psoriasis skin - GSE109248")

# ====== 2. GSE167923: CLE blood ======
print("\n\n2. CLE BLOOD - GSE167923")
print("="*60)
expr = pd.read_csv(os.path.join(DATA_DIR, 'GSE167923_expression.csv'))
expr = expr.set_index(expr.columns[0])
print(f"\nCLE blood samples: {expr.shape[1]}")
corr_matrix(expr, expr.columns.tolist(), "CLE blood - GSE167923")

# ====== 3. Compare subtype correlation strengths ======
print("\n\n3. TLR7→BAFF CORRELATION ACROSS CONDITIONS")
print("="*60)
print(f"{'Comparison':25s} {'TLR7-CXCL10':>12s} {'TLR7-BAFF':>12s} {'BAFF-IRF7':>12s} {'IRF7-CXCL10':>12s} {'TLR7-MYD88':>12s}")
print("-"*85)

# I'll need to compute these for each group
def get_pair_corr(df_genes, cols, g1, g2):
    if g1 not in df_genes.index or g2 not in df_genes.index:
        return None, None
    v1 = df_genes.loc[g1, cols].values.astype(float)
    v2 = df_genes.loc[g2, cols].values.astype(float)
    valid = ~(np.isnan(v1) | np.isnan(v2))
    if sum(valid) < 5:
        return None, None
    r, p = stats.pearsonr(v1[valid], v2[valid])
    return r, p

pairs = [('TLR7','CXCL10'), ('TLR7','TNFSF13B'), ('TNFSF13B','IRF7'), ('IRF7','CXCL10'), ('TLR7','MYD88')]
for label, samples, d in [
    ('CLE skin (109)', cle_samples, genes109),
    ('Control skin (109)', ctrl_samples, genes109),
    ('Psoriasis (109)', pso_samples, genes109),
    ('CLE blood (167)', expr.columns.tolist(), expr),
]:
    vals = []
    for g1, g2 in pairs:
        r, p = get_pair_corr(d, samples, g1, g2)
        if r is not None:
            s = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            vals.append(f"{r:+.2f}{s}")
        else:
            vals.append('N/A')
    print(f"{label:25s}  {vals[0]:>12s}  {vals[1]:>12s}  {vals[2]:>12s}  {vals[3]:>12s}  {vals[4]:>12s}")

print("\nDone!")
