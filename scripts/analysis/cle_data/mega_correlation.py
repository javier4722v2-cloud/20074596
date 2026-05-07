"""MEGA correlation analysis - all CLE + SLE data"""
import pandas as pd
import numpy as np
import os, re
from scipy import stats
from collections import defaultdict
from io import StringIO

DATA_DIR = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data'
SLE_DIR = r'C:\Users\javie\.openclaw\workspace\lumus\data'
OUT_DIR = os.path.join(DATA_DIR, 'results')
os.makedirs(OUT_DIR, exist_ok=True)

# Load probe-to-gene mapping
probe_to_gene = {}
with open(os.path.join(DATA_DIR, 'probe_to_gene.txt'), encoding='utf-8') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            probe_to_gene[parts[0]] = parts[1]

def collapse_probes(df, mapping):
    gp = defaultdict(list)
    for p in df.index:
        if p in mapping:
            gp[mapping[p]].append(p)
    ge = {}
    for g, ps in gp.items():
        ge[g] = df.loc[ps[0]] if len(ps) == 1 else df.loc[ps].median()
    return pd.DataFrame(ge).T

def parse_series_matrix(fp):
    with open(fp, encoding='utf-8', errors='replace') as f:
        t = f.read()
    lns = t.split('\n')
    titles, geo_ids = [], []
    for l in lns:
        if l.startswith('!Sample_title'): titles = [s.strip().strip('"') for s in l.split('\t')[1:]]
        if l.startswith('!Sample_geo_accession'): geo_ids = [s.strip().strip('"') for s in l.split('\t')[1:]]
    start = next(i for i,l in enumerate(lns) if 'series_matrix_table_begin' in l) + 1
    dl = []
    for i in range(start, len(lns)):
        if 'series_matrix_table_end' in lns[i]: break
        dl.append(lns[i])
    df = pd.read_csv(StringIO('\n'.join(dl)), sep='\t', index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    return df, titles, geo_ids

def corr_pair(df_genes, cols, g1, g2):
    if g1 not in df_genes.index or g2 not in df_genes.index:
        return None, None
    v1 = df_genes.loc[g1, cols].values.astype(float)
    v2 = df_genes.loc[g2, cols].values.astype(float)
    valid = ~(np.isnan(v1) | np.isnan(v2))
    if sum(valid) < 5:
        return None, None
    r, p = stats.pearsonr(v1[valid], v2[valid])
    return r, p

def analyze_group(df_g, cols, label):
    """Run correlation analysis on a group of samples"""
    pairs = [('TLR7','CXCL10'), ('TLR7','TNFSF13B'), ('TNFSF13B','IRF7'), 
             ('IRF7','CXCL10'), ('TLR7','MYD88'), ('TLR7','IRF7'),
             ('MYD88','TNFSF13B'), ('IRF7','STAT1'), ('TLR7','STAT1'),
             ('TNFSF13B','CXCL10'), ('ELANE','CXCL10'), ('IRF7','IFI44L')]
    
    results = {'label': label, 'n': len(cols)}
    for g1, g2 in pairs:
        r, p = corr_pair(df_g, cols, g1, g2)
        results[f'{g1}_{g2}'] = (r, p)
    return results

def print_group(results):
    pairs = results.keys() - {'label', 'n'}
    print(f"\n  {results['label']} (n={results['n']}):")
    for k in sorted(pairs):
        r, p = results[k]
        if r is not None:
            s = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"    {k:30s}  r={r:+.3f}  p={p:.1e} {s}")
        else:
            print(f"    {k:30s}  N/A")

# ============ COLLECT ALL DATASETS ============
all_groups = []  # (df, columns, label)

# === 1. GSE109248 - CLE skin ===
print("="*60)
print("LOADING GSE109248 (CLE skin)")
print("="*60)
df109, tits, gids = parse_series_matrix(os.path.join(DATA_DIR, 'GSE109248_series_matrix.txt'))
g109 = collapse_probes(df109, probe_to_gene)

def cls(t):
    t = t.lower()
    if 'chronic' in t: return 'CCLE'
    elif 'subacute' in t: return 'SCLE'
    elif 'acute' in t: return 'ACLE'
    elif 'psoriasis' in t: return 'Psoriasis'
    elif 'control' in t or 'normal' in t: return 'Control'
    return 'Other'

grps109 = [cls(t) for t in tits]
sg109 = defaultdict(list)
for gid, g in zip(gids, grps109): sg109[g].append(gid)

for lbl, cols in [
    ('CLE skin (GSE109248)', sg109['CCLE']+sg109['SCLE']+sg109['ACLE']),
    ('CLE-CCLE skin', sg109['CCLE']),
    ('CLE-SCLE skin', sg109['SCLE']),
    ('CLE-ACLE skin', sg109['ACLE']),
    ('Control skin', sg109['Control']),
    ('Psoriasis skin', sg109['Psoriasis']),
]:
    all_groups.append((g109, cols, lbl))

print(f"  CLE: {len(sg109['CCLE']+sg109['SCLE']+sg109['ACLE'])}")
print(f"  Control: {len(sg109['Control'])}")
print(f"  Psoriasis: {len(sg109['Psoriasis'])}")

# === 2. GSE112943 - CLE skin replicate ===
print("\n" + "="*60)
print("LOADING GSE112943 (CLE skin replicate)")
print("="*60)
df112, tits112, gids112 = parse_series_matrix(os.path.join(DATA_DIR, 'GSE112943_series_matrix.txt'))
g112 = collapse_probes(df112, probe_to_gene)

grps112 = [cls(t) for t in tits112]
sg112 = defaultdict(list)
for gid, g in zip(gids112, grps112): sg112[g].append(gid)

for lbl, cols in [
    ('CLE skin (GSE112943)', sg112['CCLE']+sg112['SCLE']),
    ('CLE-CCLE skin rep', sg112['CCLE']),
    ('CLE-SCLE skin rep', sg112['SCLE']),
    ('Control skin rep', sg112['Control']),
]:
    all_groups.append((g112, cols, lbl))

print(f"  CLE: {len(sg112['CCLE']+sg112['SCLE'])}")
print(f"  Control: {len(sg112['Control'])}")
# LN has no control kidney, skip
print(f"  LN kidney (no control): {len(sg112['LN'])}")

# === 3. GSE167923 - CLE blood ===
print("\n" + "="*60)
print("LOADING GSE167923 (CLE blood)")
print("="*60)
expr167 = pd.read_csv(os.path.join(DATA_DIR, 'GSE167923_expression.csv'))
expr167 = expr167.set_index(expr167.columns[0])
all_groups.append((expr167, list(expr167.columns), 'CLE blood (GSE167923)'))
print(f"  CLE blood: {expr167.shape[1]}")

# === 4. GSE81622 - SLE blood (Illumina GPL10558) ===
print("\n" + "="*60)
print("LOADING GSE81622 (SLE blood, GPL10558)")
print("="*60)
# Parse non-normalized.txt to extract AVG_Signal columns
sle_file = os.path.join(SLE_DIR, 'gse81622', 'GSE81622_non-normalized.txt')
# Read only AVG_Signal columns
with open(sle_file, encoding='utf-8', errors='replace') as f:
    header = f.readline().strip().split('\t')

# Find AVG_Signal columns
avg_cols = [i for i, h in enumerate(header) if 'AVG_Signal' in h]
probe_col = 0  # ID_REF
print(f"  Found {len(avg_cols)} AVG_Signal columns out of {len(header)} total")

# Read only probe ID and AVG_Signal columns
df816 = pd.read_csv(sle_file, sep='\t', usecols=[probe_col] + avg_cols, index_col=0,
                    encoding='utf-8', low_memory=False)
# Clean column names to just sample IDs
new_cols = []
for c in df816.columns:
    # Extract sample prefix: e.g. "9461595024_A.AVG_Signal" -> keep just the column
    new_cols.append(c.replace('.AVG_Signal', ''))
df816.columns = new_cols

# Collapse to genes
g816 = collapse_probes(df816, probe_to_gene)
print(f"  Probes: {df816.shape[0]} → Genes: {g816.shape[0]}, Samples: {g816.shape[1]}")

# Load metadata to split SLE vs Control
meta_path = os.path.join(SLE_DIR, 'gse81622', 'GSE81622_metadata.txt')
if os.path.exists(meta_path):
    meta = pd.read_csv(meta_path, sep='\t')
    # Determine condition from metadata
    sle_cols = []
    ctrl_cols = []
    for col in g816.columns:
        # Try to match to metadata
        pass
    
    # Try a different approach: classify by sample name prefix convention
    # GSE81622: first group = healthy controls, second group = SLE
    # From the report, we know the groups
    print(f"  Metadata found: {meta.shape}")
    print(f"  First 5 samples: {list(g816.columns[:5])}")

# For now, split by scanning sample IDs against metadata
meta_text = open(meta_path, encoding='utf-8').read()
sle_samples = []
ctrl_samples = []
for col in g816.columns:
    gsm_id = col.split('_')[0] if '_' in col else col
    if gsm_id in meta_text:
        # Find the line
        for line in meta_text.split('\n'):
            if gsm_id in line and 'SLE' in line:
                sle_samples.append(col)
                break
            elif gsm_id in line and 'Control' in line:
                ctrl_samples.append(col)
                break

if sle_samples or ctrl_samples:
    # We found some matches
    all_groups.append((g816, sle_samples, 'SLE blood (GSE81622)'))
    all_groups.append((g816, ctrl_samples, 'Control blood (GSE81622)'))
    print(f"  SLE blood: {len(sle_samples)}, Control blood: {len(ctrl_samples)}")
else:
    # Alternative: split by index (first half = control, second = SLE based on typical GEO format)
    n = len(g816.columns)
    half = n // 2
    all_groups.append((g816, list(g816.columns[half:]), 'SLE blood (GSE81622)'))
    all_groups.append((g816, list(g816.columns[:half]), 'Control blood (GSE81622)'))
    print(f"  Auto-split: SLE {n-half}, Control {half}")

# === 5. GSE50772 - SLE blood (Affymetrix) ===
print("\n" + "="*60)
print("LOADING GSE50772 (SLE blood, Affymetrix)")
print("="*60)
expr507 = pd.read_csv(os.path.join(os.path.dirname(DATA_DIR), 'GSE50772_expression_matrix.csv'))
probe_col = expr507.columns[0]
expr507 = expr507.set_index(probe_col)

# Load sample info
si = pd.read_csv(os.path.join(os.path.dirname(DATA_DIR), 'GSE50772_sample_info.csv'), index_col=0)
sle507 = [c for c in expr507.columns if c in si.index and si.loc[c, 'condition'] == 'SLE']
ctrl507 = [c for c in expr507.columns if c in si.index and si.loc[c, 'condition'] == 'Control']
print(f"  SLE blood: {len(sle507)}, Control: {len(ctrl507)}")
all_groups.append((expr507, sle507, 'SLE blood (GSE50772)'))
all_groups.append((expr507, ctrl507, 'Control blood (GSE50772)'))

# ============ RUN CORRELATION ON ALL GROUPS ============
print("\n\n" + "="*70)
print("CORRELATION ANALYSIS - ALL GROUPS")
print("="*70)

results_list = []
for df_g, cols, label in all_groups:
    if len(cols) < 3:
        continue
    res = analyze_group(df_g, cols, label)
    results_list.append(res)

# Print key pairs across all groups
print("\n\n" + "="*70)
print("SUMMARY: TLR7→CXCL10 CORRELATION ACROSS ALL CONDITIONS")
print("="*70)
for res in sorted(results_list, key=lambda r: r.get('n', 0), reverse=True):
    r, p = res.get('TLR7_CXCL10', (None, None))
    s = ''
    if r is not None:
        s = f"{r:+.3f}" + ('*'*3 if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else '')
    print(f"  {res['label']:35s}  n={res['n']:3d}  TLR7-CXCL10: {s}")

print("\n\nSUMMARY: TLR7→BAFF CORRELATION")
for res in sorted(results_list, key=lambda r: r.get('n', 0), reverse=True):
    r, p = res.get('TLR7_TNFSF13B', (None, None))
    s = ''
    if r is not None:
        s = f"{r:+.3f}" + ('*'*3 if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else '')
    print(f"  {res['label']:35s}  n={res['n']:3d}  TLR7-BAFF:   {s}")

print("\n\nSUMMARY: BAFF→IRF7 CORRELATION")
for res in sorted(results_list, key=lambda r: r.get('n', 0), reverse=True):
    r, p = res.get('TNFSF13B_IRF7', (None, None))
    s = ''
    if r is not None:
        s = f"{r:+.3f}" + ('*'*3 if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else '')
    print(f"  {res['label']:35s}  n={res['n']:3d}  BAFF-IRF7:   {s}")

print("\n\nSUMMARY: IRF7→CXCL10 CORRELATION")
for res in sorted(results_list, key=lambda r: r.get('n', 0), reverse=True):
    r, p = res.get('IRF7_CXCL10', (None, None))
    s = ''
    if r is not None:
        s = f"{r:+.3f}" + ('*'*3 if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else '')
    print(f"  {res['label']:35s}  n={res['n']:3d}  IRF7-CXCL10: {s}")

print("\n\nSUMMARY: ELANE→CXCL10 (NETs signature)")
for res in sorted(results_list, key=lambda r: r.get('n', 0), reverse=True):
    r, p = res.get('ELANE_CXCL10', (None, None))
    s = ''
    if r is not None:
        s = f"{r:+.3f}" + ('*'*3 if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else '')
    print(f"  {res['label']:35s}  n={res['n']:3d}  ELANE-CXCL10:{s}")

# Save to CSV
summary_rows = []
headers = ['label', 'n', 'TLR7-CXCL10', 'TLR7-BAFF', 'BAFF-IRF7', 'IRF7-CXCL10', 'TLR7-MYD88', 'ELANE-CXCL10', 'IRF7-IFI44L']
print("\n\n" + "="*70)
print("FULL TABLE")
print("="*70)
print(f"{'Group':35s} {'n':>4s} {'TLR7→CXCL10':>12s} {'TLR7→BAFF':>12s} {'BAFF→IRF7':>12s} {'IRF7→CXCL10':>12s} {'TLR7→MYD88':>12s} {'ELANE→CXCL10':>12s} {'IRF7→IFI44L':>12s}")
print("-"*125)

pairs_for_table = [('TLR7_CXCL10', 'TLR7→CXCL10'), ('TLR7_TNFSF13B', 'TLR7→BAFF'),
                   ('TNFSF13B_IRF7', 'BAFF→IRF7'), ('IRF7_CXCL10', 'IRF7→CXCL10'),
                   ('TLR7_MYD88', 'TLR7→MYD88'), ('ELANE_CXCL10', 'ELANE→CXCL10'),
                   ('IRF7_IFI44L', 'IRF7→IFI44L')]

# Sort by category
def group_order(lbl):
    if 'SLE blood' in lbl: return 0
    if 'CLE blood' in lbl: return 1
    if 'CLE skin' in lbl: return 2
    if 'Psoriasis' in lbl: return 4
    if 'Control skin' in lbl: return 5
    if 'Control blood' in lbl: return 6
    if 'CCLE' in lbl: return 3
    if 'SCLE' in lbl: return 3
    if 'ACLE' in lbl: return 3
    return 9
results_list.sort(key=lambda r: (group_order(r['label']), -r.get('n', 0)))

for res in results_list:
    vals = [res['label'], res['n']]
    for key, _ in pairs_for_table:
        r, p = res.get(key, (None, None))
        if r is not None:
            s = '*'*3 if p<0.001 else '*'*2 if p<0.01 else '*' if p<0.05 else ''
            vals.append(f"{r:+.2f}{s}")
        else:
            vals.append('N/A')
    
    print(f"{vals[0]:35s} {vals[1]:4d} {vals[2]:>12s} {vals[3]:>12s} {vals[4]:>12s} {vals[5]:>12s} {vals[6]:>12s} {vals[7]:>12s} {vals[8]:>12s}")

print("\nDone! All groups analyzed.")
