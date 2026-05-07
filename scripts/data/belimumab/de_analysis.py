"""
Belimumab SLE Dataset: Paired Differential Expression Analysis
PRE (0M) vs POST (6M) treatment, paired t-test per gene
"""
import pandas as pd
import numpy as np
from scipy import stats
import re
import json
import time
import warnings
warnings.filterwarnings('ignore')

DATA_PATH = r"C:\Users\javie\.openclaw\workspace\lumus\data\belimumab\SLE_BLM_normalized_count.txt"
OUT_CSV = r"C:\Users\javie\.openclaw\workspace\lumus\data\belimumab\full_de.csv"
OUT_SIG = r"C:\Users\javie\.openclaw\workspace\lumus\data\belimumab\significant_genes.json"

print("=" * 70)
print("BELLMUMAB SLE: Paired Differential Expression Analysis")
print("=" * 70)
print()

# ── 1. Load data ──────────────────────────────────────────────────────────────
print("📂 Loading data...")
t0 = time.time()
df = pd.read_csv(DATA_PATH, sep='\t', low_memory=False)
print(f"   Loaded {df.shape[0]} genes × {df.shape[1]} columns ({df.shape[1]-1} samples)")
print(f"   Time: {time.time()-t0:.1f}s")
print()

# ── 2. Parse columns ──────────────────────────────────────────────────────────
cols = df.columns.tolist()
print("📋 Column categories:")
for c in cols[1:]:
    # parse sample name: e.g. "R_0M_10901"
    parts = c.split('_')
    # Standardize: original had "CR_" and "CNR_" prefixes, but clean split shows "R_" and "NR_"
    pass  # handled below

# Count by type
healthy_cols = [c for c in cols[1:] if c.startswith('healthy')]
r_cols = [c for c in cols[1:] if c.startswith('R_') or c.lower().startswith('cr_')]
nr_cols = [c for c in cols[1:] if c.startswith('NR_') or c.lower().startswith('cnr_')]

# Renormalize: some columns might have "CR_" / "CNR_" prefix in the raw file
all_sle_cols = r_cols + nr_cols

print(f"   Healthy controls: {len(healthy_cols)}")
print(f"   Responders (R):   {len(r_cols)}")
print(f"   Non-Responders:   {len(nr_cols)}")
print(f"   Total SLE:        {len(all_sle_cols)}")
print()

# ── 3. Parse timepoints ──────────────────────────────────────────────────────
def parse_sample(s):
    """Return (group, timepoint, patient_id) from sample name."""
    s_clean = s.replace('CR_', 'R_').replace('CNR_', 'NR_')
    # Handle healthy separately
    if s_clean.startswith('healthy'):
        match = re.match(r'healthy_(\d+)', s_clean)
        if match:
            return ('healthy', '0M', match.group(1))
        return ('healthy', '0M', s_clean)
    
    match = re.match(r'(R|NR)_(\d+M)_(\d+)', s_clean)
    if match:
        return (match.group(1), match.group(2), match.group(3))
    return None

parsed = {}
for c in cols[1:]:
    p = parse_sample(c)
    if p:
        parsed[c] = p
    else:
        print(f"   ⚠️  Could not parse: {c}")

# ── 4. Build patient pairs (PRE=0M, POST=6M or 3M) ───────────────────────────
# Map patient base ID to samples at each timepoint
patients = {}  # {base_id: {group: 'R'|'NR', 0M: col, 3M: col, 6M: col}}

for c, (grp, tp, pid) in parsed.items():
    if grp == 'healthy':
        continue
    # Extract base patient ID (remove trailing visit digit)
    # Pattern: pid could be like 10901 (base=1090, visit=1)
    base_match = re.match(r'(\d+)(\d)$', pid)
    if base_match:
        base_id = base_match.group(1)
    else:
        base_id = pid
    
    key = f"{grp}_{base_id}"
    if key not in patients:
        patients[key] = {'group': grp, 'base_id': base_id}
    patients[key][tp] = c

# Count paired patients
paired_0m_6m = {k: v for k, v in patients.items() if '0M' in v and '6M' in v}
paired_0m_3m = {k: v for k, v in patients.items() if '0M' in v and '3M' in v and '6M' not in v}
print(f"🔗 Paired patients (0M + 6M): {len(paired_0m_6m)}")
print(f"🔗 Paired patients (0M + 3M): {len(paired_0m_3m)}")

r_paired = {k: v for k, v in paired_0m_6m.items() if v['group'] == 'R'}
nr_paired = {k: v for k, v in paired_0m_6m.items() if v['group'] == 'NR'}
print(f"   Responders paired: {len(r_paired)}")
print(f"   Non-Responders paired: {len(nr_paired)}")
print()

# ── 5. Prepare paired data matrices ──────────────────────────────────────────
# Set gene as index
df_indexed = df.set_index('Gene')

# Build PRE and POST matrices
pre_cols = [v['0M'] for v in paired_0m_6m.values()]
post_cols = [v['6M'] for v in paired_0m_6m.values()]

# Add 3M patients that don't have 6M
for k, v in paired_0m_3m.items():
    pre_cols.append(v['0M'])
    post_cols.append(v['3M'])

all_paired = {**paired_0m_6m, **paired_0m_3m}
pre_cols_all = [v['0M'] for v in all_paired.values()]
post_cols_all = [v[t] for v, t in zip(all_paired.values(), 
                                       ['6M' if '6M' in v else '3M' for v in all_paired.values()])]

print(f"📊 Total paired samples for DE: {len(pre_cols_all)}")
print(f"   PRE cols: {pre_cols_all[:3]}...{pre_cols_all[-1]}")
print(f"   POST cols: {post_cols_all[:3]}...{post_cols_all[-1]}")
print()

# Extract expression matrices
pre_data = df_indexed[pre_cols_all].values.astype(np.float64)
post_data = df_indexed[post_cols_all].values.astype(np.float64)
gene_names = df_indexed.index.values

# Filter out genes that are all zero in both conditions
non_zero_mask = ~((pre_data.sum(axis=1) == 0) & (post_data.sum(axis=1) == 0))
pre_data = pre_data[non_zero_mask]
post_data = post_data[non_zero_mask]
gene_names = gene_names[non_zero_mask]
n_genes = len(gene_names)
print(f"🧬 Genes after removing all-zero: {n_genes}")
print()

# ── 6. Paired t-test ─────────────────────────────────────────────────────────
print("⚗️  Running paired t-tests...")
t0 = time.time()

results = []
n_pairs = pre_data.shape[1]

for i in range(n_genes):
    pre = pre_data[i, :]
    post = post_data[i, :]
    
    # Only use pairs where both values are non-zero
    valid = (pre > 0) & (post > 0)
    n_valid = valid.sum()
    
    if n_valid >= 3:  # Need at least 3 pairs for a meaningful test
        pre_v = pre[valid]
        post_v = post[valid]
        
        t_stat, p_val = stats.ttest_rel(pre_v, post_v)
        
        # Log2 fold change: mean(post) - mean(pre)  (positive = up after treatment)
        mean_pre = np.mean(pre_v)
        mean_post = np.mean(post_v)
        if mean_pre > 0 and mean_post > 0:
            log2fc = np.log2(mean_post) - np.log2(mean_pre)
        elif mean_post > 0 and mean_pre == 0:
            log2fc = np.log2(mean_post)  # from near-zero
        elif mean_pre > 0 and mean_post == 0:
            log2fc = -np.log2(mean_pre)  # down to zero
        else:
            log2fc = 0
    else:
        t_stat = 0
        p_val = 1.0
        log2fc = 0
        n_valid = 0
    
    results.append({
        'Gene': gene_names[i],
        'log2FC': log2fc,
        't_statistic': t_stat,
        'p_value': p_val,
        'mean_pre': mean_pre if 'mean_pre' in dir() else 0,
        'mean_post': mean_post if 'mean_post' in dir() else 0,
        'n_pairs_valid': n_valid
    })

results_df = pd.DataFrame(results)

# Multiple testing correction (Bonferroni)
results_df['p_adjusted_bonf'] = np.minimum(results_df['p_value'] * n_genes, 1.0)

# FDR (Benjamini-Hochberg)
from scipy.stats import rankdata
p_vals = results_df['p_value'].values
ranked = rankdata(p_vals)
n_tests = len(p_vals)
results_df['p_adjusted_fdr'] = np.minimum(p_vals * n_tests / ranked, 1.0)

# Sort by p-value
results_df = results_df.sort_values('p_value').reset_index(drop=True)

print(f"   Done in {time.time()-t0:.1f}s")
print()

# ── 7. Save full results ─────────────────────────────────────────────────────
results_df.to_csv(OUT_CSV, index=False)
print(f"💾 Full DE results saved: {OUT_CSV}")

# ── 8. Significant genes ──────────────────────────────────────────────────────
sig = results_df[results_df['p_value'] < 0.05].copy()
sig_genes = sig['Gene'].tolist()
print(f"📈 Significant genes (p<0.05): {len(sig_genes)} / {n_genes} ({len(sig_genes)/n_genes*100:.1f}%)")

# FDR significant
sig_fdr = results_df[results_df['p_adjusted_fdr'] < 0.05].copy()
print(f"📈 Significant (FDR<0.05): {len(sig_fdr)} / {n_genes}")

sig_json = {
    'metadata': {
        'total_genes_analyzed': n_genes,
        'significant_p_0.05': len(sig_genes),
        'significant_fdr_0.05': len(sig_fdr),
        'n_paired_samples': n_pairs,
        'n_responders_paired': len(r_paired),
        'n_nonresponders_paired': len(nr_paired),
        'analysis_type': 'paired_t_test_pre_0m_vs_post_6m_3m'
    },
    'significant_genes': sig_genes,
    'significant_fdr_genes': sig_fdr['Gene'].tolist() if len(sig_fdr) > 0 else []
}

with open(OUT_SIG, 'w') as f:
    json.dump(sig_json, f, indent=2)
print(f"💾 Significant genes saved: {OUT_SIG}")
print()

# ── 9. Print summary ─────────────────────────────────────────────────────────
print("=" * 70)
print("📊 RESULTS SUMMARY")
print("=" * 70)
print()

# Top 20 upregulated (post > pre, log2FC > 0)
top_up = results_df[results_df['log2FC'] > 0].head(20)
print("🔺 TOP 20 UPREGULATED after Belimumab:")
print(f"{'Gene':<35} {'log2FC':>8} {'p_value':>10} {'p_adj(FDR)':>10}")
print("-" * 65)
for _, r in top_up.iterrows():
    gene_short = r['Gene'][:34] if len(str(r['Gene'])) > 34 else r['Gene']
    print(f"{gene_short:<35} {r['log2FC']:>8.4f} {r['p_value']:>10.2e} {r['p_adjusted_fdr']:>10.2e}")

print()

# Top 20 downregulated (post < pre, log2FC < 0)
top_down = results_df[results_df['log2FC'] < 0].head(20)
print("🔻 TOP 20 DOWNREGULATED after Belimumab:")
print(f"{'Gene':<35} {'log2FC':>8} {'p_value':>10} {'p_adj(FDR)':>10}")
print("-" * 65)
for _, r in top_down.iterrows():
    gene_short = r['Gene'][:34] if len(str(r['Gene'])) > 34 else r['Gene']
    print(f"{gene_short:<35} {r['log2FC']:>8.4f} {r['p_value']:>10.2e} {r['p_adjusted_fdr']:>10.2e}")

print()

# ── 10. Functional categorization ─────────────────────────────────────────────
print("=" * 70)
print("🧬 FUNCTIONAL CATEGORIZATION")
print("=" * 70)
print()

# Interferon signature genes
ifn_genes_list = ['IFIT1', 'IFI44', 'MX1', 'OAS', 'ISG15', 'IFIT2', 'IFIT3', 
                  'IFI44L', 'IFITM1', 'IFITM3', 'MX2', 'OAS1', 'OAS2', 'OAS3',
                  'OASL', 'IRF7', 'IRF9', 'STAT1', 'STAT2', 'USP18', 'RSAD2',
                  'ISG20', 'HERC5', 'HERC6', 'BST2', 'IFI16', 'IFIH1', 'IFIT5']
# BAFF / B-cell related
baff_genes = ['TNFSF13B', 'TNFRSF13B', 'TNFRSF13C', 'MS4A1', 'CD19', 'CD22',
              'CD79A', 'CD79B', 'PAX5', 'BLK', 'BACH2', 'BCL6', 'BANK1']
# TLR / MYD88
tlr_genes = ['TLR1', 'TLR2', 'TLR3', 'TLR4', 'TLR5', 'TLR6', 'TLR7', 'TLR8',
             'TLR9', 'TLR10', 'MYD88', 'IRAK1', 'IRAK4', 'TRAF6']
# JAK/STAT
jak_stat_genes = ['JAK1', 'JAK2', 'JAK3', 'TYK2', 'STAT1', 'STAT2', 'STAT3',
                  'STAT4', 'STAT5A', 'STAT5B', 'STAT6']
# Cytokines
cyto_genes = ['IL1B', 'IL2', 'IL4', 'IL5', 'IL6', 'IL10', 'IL12A', 'IL12B',
              'IL17A', 'IL17F', 'IL21', 'IL23A', 'TNF', 'TNFRSF1A',
              'IFNG', 'TGFB1', 'CCL2', 'CCL5', 'CXCL10', 'CXCL9', 'CXCL8',
              'CX3CL1', 'CCL20', 'CCL19', 'CCL21', 'CCL17', 'CCL22']

categories = {
    'Interferon_signature': ifn_genes_list,
    'BAFF_Bcell': baff_genes,
    'TLR_MYD88': tlr_genes,
    'JAK_STAT': jak_stat_genes,
    'Cytokines': cyto_genes
}

# Create a gene lookup: extract gene symbol from ENSG..._SYMBOL
def extract_symbol(gene_id):
    if '_' in gene_id:
        return gene_id.split('_', 1)[1] if gene_id.startswith('ENSG') else gene_id
    return gene_id

results_df['symbol'] = results_df['Gene'].apply(extract_symbol)

for cat_name, gene_list in categories.items():
    print(f"── {cat_name} ──")
    cat_results = results_df[results_df['symbol'].str.upper().isin([g.upper() for g in gene_list])].copy()
    
    if len(cat_results) == 0:
        print("   (no genes found in dataset)")
        print()
        continue
    
    cat_sig = cat_results[cat_results['p_value'] < 0.05]
    n_sig = len(cat_sig)
    n_total = len(cat_results)
    
    print(f"   Found {n_total} genes, {n_sig} significant (p<0.05)")
    
    if n_total > 0:
        avg_fc = cat_results['log2FC'].mean()
        direction = "UP" if avg_fc > 0 else "DOWN"
        print(f"   Average log2FC: {avg_fc:.4f} ({direction}-regulated overall)")
    
    if n_sig > 0:
        print(f"\n   {'Gene':<35} {'log2FC':>8} {'p_value':>10} {'Direction':>10}")
        print("   " + "-" * 65)
        for _, r in cat_sig.sort_values('p_value').iterrows():
            dir_str = "UP ↑" if r['log2FC'] > 0 else "DOWN ↓"
            gene_short = r['symbol'][:34] if len(str(r['symbol'])) > 34 else r['symbol']
            print(f"   {gene_short:<35} {r['log2FC']:>8.4f} {r['p_value']:>10.2e} {dir_str:>10}")
    print()

# ── 11. Key specific questions ────────────────────────────────────────────────
print("=" * 70)
print("🔑 KEY QUESTIONS")
print("=" * 70)
print()

# BAFF-R (TNFRSF13C)
baffr = results_df[results_df['symbol'].str.upper() == 'TNFRSF13C']
if len(baffr) > 0:
    r = baffr.iloc[0]
    dir_str = "DOWN (disminuye)" if r['log2FC'] < 0 else "UP (aumenta)"
    sig_str = "SIGNIFICATIVO" if r['p_value'] < 0.05 else "no significativo"
    print(f"❓ ¿BAFF-R (TNFRSF13C) baja después de Belimumab?")
    print(f"   → log2FC = {r['log2FC']:.4f} ({dir_str})")
    print(f"   → p-value = {r['p_value']:.2e} ({sig_str})")
else:
    print(f"❓ ¿BAFF-R (TNFRSF13C)?")
    print("   → No encontrado en el dataset")
print()

# Interferon signature summary
ifn_sig_cat = results_df[results_df['symbol'].str.upper().isin([g.upper() for g in ifn_genes_list])]
if len(ifn_sig_cat) > 0:
    avg_fc_ifn = ifn_sig_cat['log2FC'].mean()
    p_vals_ifn = ifn_sig_cat['p_value'].dropna()
    n_ifn_sig = (ifn_sig_cat['p_value'] < 0.05).sum()
    n_ifn_total = len(ifn_sig_cat)
    
    fisher_method_stat = -2 * np.sum(np.log(np.clip(p_vals_ifn, 1e-300, None)))
    combined_p = stats.chi2.sf(fisher_method_stat, 2 * len(p_vals_ifn))
    
    dir_str = "disminuye" if avg_fc_ifn < 0 else "aumenta"
    print(f"❓ ¿La firma interferón disminuye después de Belimumab?")
    print(f"   → {n_ifn_sig}/{n_ifn_total} genes IFN significativos")
    print(f"   → Media log2FC de firma IFN: {avg_fc_ifn:.4f} (globalmente {dir_str})")
    print(f"   → Combined p-value (Fisher): {combined_p:.2e}")
    
    if avg_fc_ifn < 0 and combined_p < 0.05:
        print(f"   ✅ CONCLUSIÓN: La firma interferón DISMINUYE significativamente")
        print(f"      después del tratamiento con Belimumab")
    elif avg_fc_ifn < 0:
        print(f"   ⚠️  Tendencia a disminución pero no significativa (Fisher p={combined_p:.3f})")
    else:
        print(f"   ⚠️  La firma interferón no muestra disminución clara")
else:
    print("❓ ¿La firma interferón disminuye?")
    print("   → Genes IFN no encontrados en el dataset")
print()

print("=" * 70)
print("✅ ANALYSIS COMPLETE")
print("=" * 70)
