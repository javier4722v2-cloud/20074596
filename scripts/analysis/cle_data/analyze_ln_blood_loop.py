import gzip
import numpy as np
from scipy import stats
import json

# Load GSE72747 (LN blood, 30 samples)
f = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data\GSE72747_series_matrix.txt.gz'
c = gzip.open(f, 'rt', encoding='utf-8').read()
lines = c.split('\n')

# Parse expression data
expr = {}
header = None
in_table = False
for l in lines:
    if '!series_matrix_table_begin' in l:
        in_table = True
        continue
    if '!series_matrix_table_end' in l:
        break
    if in_table:
        if header is None:
            header = l.strip().split('\t')
        else:
            parts = l.strip().split('\t')
            probe = parts[0].strip('"')
            values = np.array([float(v.strip('"')) for v in parts[1:]])
            expr[probe] = values

# Target genes with their probes
gene_probes = {
    'TLR7': ['220146_at'],
    'MYD88': ['209124_at'],
    'BAFF': ['223500_at'],  # TNFSF13B
    'IRF7': ['208436_s_at'],
    'CXCL10': ['204533_at'],
    'ELANE': ['206871_at'],
    'IFI27': ['202411_at'],
    'IFI44L': ['214453_s_at'],
    'IFIT1': ['203153_at'],
    'ISG15': ['205483_s_at'],
    'MX1': ['202086_at'],
    'TNF': ['207113_s_at'],
    'IL6': ['205207_at'],
    'STAT1': ['209969_s_at'],
}

# Get gene-level expression (single probe per gene)
gene_expr = {}
for gene, probes in gene_probes.items():
    probe = probes[0]
    if probe in expr:
        gene_expr[gene] = expr[probe]
        print(f'{gene}: {probe} mean={np.mean(expr[probe]):.4f} std={np.std(expr[probe]):.4f}')
    else:
        print(f'{gene}: {probe} NOT FOUND')

print(f'\nComputing correlations in LN blood (n=30)...')
print(f'{"Pair":<25} {"r":<8} {"p":<10} {"Sig":<6}')
print('-' * 50)

loop_pairs = [
    ('TLR7', 'MYD88'),
    ('MYD88', 'IRF7'),
    ('TLR7', 'BAFF'),
    ('BAFF', 'IRF7'),
    ('IRF7', 'CXCL10'),
    ('IRF7', 'IFIT1'),
    ('IRF7', 'ISG15'),
    ('TLR7', 'CXCL10'),
    ('ELANE', 'CXCL10'),
    ('TLR7', 'ELANE'),
    ('ELANE', 'TLR7'),
]

results = {}
for g1, g2 in loop_pairs:
    if g1 in gene_expr and g2 in gene_expr:
        r, p = stats.pearsonr(gene_expr[g1], gene_expr[g2])
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        print(f'{g1} vs {g2:<12} {r:>+.4f}  {p:<8.6f}  {sig}')
        results[f'{g1}_vs_{g2}'] = {'r': round(r, 4), 'p': round(p, 6), 'sig': sig}

# Save results
out = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data\GSE72747_LN_blood_correlations.json'
with open(out, 'w') as fh:
    json.dump(results, fh, indent=2)
print(f'\nSaved to {out}')
