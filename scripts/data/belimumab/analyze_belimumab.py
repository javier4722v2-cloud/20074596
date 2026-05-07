import pandas as pd
import numpy as np
import json

f = r'C:\Users\javie\.openclaw\workspace\lumus\data\belimumab\SLE_BLM_normalized_count.txt'
out = r'C:\Users\javie\.openclaw\workspace\lumus\belimumab_validation.md'

# Load
print("Loading data...")
df = pd.read_csv(f, sep='\t')
gene_col = 'Gene'
print(f"Shape: {df.shape}")

# Parse genes - extract symbol after underscore
def extract_symbol(g):
    parts = g.split('_', 1)
    return parts[1] if len(parts) > 1 else parts[0]

df['Symbol'] = df[gene_col].apply(extract_symbol)

# Define target genes (exact symbol match)
target_genes = {
    'TLR7': None, 'MYD88': None, 'TNFSF13B': None, 'IRF7': None,
    'TNF': None, 'ESR1': None, 'IL6': None,
    'IFIT1': None, 'CXCL10': None, 'MX1': None,
    'IFIH1': None, 'OAS1': None, 'STAT1': None, 'STAT3': None,
    'IL6R': None, 'IL6ST': None, 'TNFRSF13B': None, 'TNFRSF13C': None,
    'TNFRSF1A': None, 'TNFSF13': None
}

# Map samples to groups
columns = [c for c in df.columns if c not in [gene_col, 'Symbol']]
groups = {}
for c in columns:
    if c.startswith('healthy_'):
        groups[c] = 'HC'
    elif '_0M_' in c:
        groups[c] = 'PRE_BELIMUMAB'
    elif '_3M_' in c:
        groups[c] = '3M_POST'
    elif '_6M_' in c:
        groups[c] = '6M_POST'
    else:
        groups[c] = 'UNKNOWN'

print(f"\nGroup sizes:")
for g in sorted(set(groups.values())):
    cnt = sum(1 for v in groups.values() if v == g)
    print(f"  {g}: {cnt}")

# Compute mean per group
group_cols = {}
for g in sorted(set(groups.values())):
    cols = [c for c in columns if groups[c] == g]
    group_cols[g] = cols

results = []
for g in target_genes:
    matches = df[df['Symbol'] == g]
    if len(matches) == 0:
        print(f"  {g}: NOT FOUND")
        continue
    if len(matches) > 1:
        # For duplicated symbols, keep first
        pass
    row = matches.iloc[0]
    gene_full = row[gene_col]
    
    entry = {'Gene': g, 'Full_ID': gene_full}
    for grp, cols in group_cols.items():
        vals = row[cols].values.astype(float)
        entry[f'mean_{grp}'] = np.mean(vals)
        entry[f'std_{grp}'] = np.std(vals)
    
    results.append(entry)

res_df = pd.DataFrame(results)

# Calculate fold changes
print("\n=== FOLD CHANGES ===")
print(f"{'Gene':<15} {'SLEvsHC':<12} {'3MvsPRE':<12} {'6MvsPRE':<12} {'6MvsHC':<12} {'Direction':<20}")
print("-"*75)

fc_data = []
for _, r in res_df.iterrows():
    sle = r['mean_PRE_BELIMUMAB']
    hc = r['mean_HC']
    m3 = r['mean_3M_POST']
    m6 = r['mean_6M_POST']
    
    fc_sle_hc = sle / hc if hc > 0 else np.nan
    fc_m3_pre = m3 / sle if sle > 0 else np.nan
    fc_m6_pre = m6 / sle if sle > 0 else np.nan
    fc_m6_hc = m6 / hc if hc > 0 else np.nan
    
    # Determine direction of change
    if fc_sle_hc > 1.15:
        sle_dir = "UP in SLE"
    elif fc_sle_hc < 0.87:
        sle_dir = "DOWN in SLE"
    else:
        sle_dir = "flat"
    
    if fc_m6_pre < 0.87:
        bel_dir = "BELIMUMAB DOWN"
    elif fc_m6_pre > 1.15:
        bel_dir = "BELIMUMAB UP"
    else:
        bel_dir = "unchanged"
    
    direction = f"{sle_dir} → {bel_dir}"
    
    print(f"{r['Gene']:<15} {fc_sle_hc:<12.3f} {fc_m3_pre:<12.3f} {fc_m6_pre:<12.3f} {fc_m6_hc:<12.3f} {direction:<20}")
    
    fc_data.append({
        'Gene': r['Gene'],
        'SLE_vs_HC': round(fc_sle_hc, 3),
        '3M_vs_PRE': round(fc_m3_pre, 3),
        '6M_vs_PRE': round(fc_m6_pre, 3),
        '6M_vs_HC': round(fc_m6_hc, 3),
        'Direction': direction,
        'Mean_HC': round(r['mean_HC'], 2),
        'Mean_PRE_BELIMUMAB': round(r['mean_PRE_BELIMUMAB'], 2),
        'Mean_3M_POST': round(r['mean_3M_POST'], 2),
        'Mean_6M_POST': round(r['mean_6M_POST'], 2)
    })

# ===== PAIRED ANALYSIS: Only patients with data at all 3 timepoints =====
print("\n\n=== PAIRED ANALYSIS (patients with PRE, 3M, and 6M) ===")
# Find patient IDs that appear in all 3 timepoints
pre_ids = set(c.split('_')[2] for c in columns if '_0M_' in c)
m3_ids = set(c.split('_')[2] for c in columns if '_3M_' in c)
m6_ids = set(c.split('_')[2] for c in columns if '_6M_' in c)
common_ids = pre_ids & m3_ids & m6_ids
print(f"Patients with all 3 timepoints: {len(common_ids)}")

paired_fc = []
for g in target_genes:
    matches = df[df['Symbol'] == g]
    if len(matches) == 0:
        continue
    row = matches.iloc[0]
    
    pre_vals = []
    m3_vals = []
    m6_vals = []
    
    for pid in common_ids:
        pre_col = [c for c in columns if f'_{pid}' in c and '_0M_' in c]
        m3_col = [c for c in columns if f'_{pid}' in c and '_3M_' in c]
        m6_col = [c for c in columns if f'_{pid}' in c and '_6M_' in c]
        
        if pre_col and m3_col and m6_col:
            pre_vals.append(float(row[pre_col[0]]))
            m3_vals.append(float(row[m3_col[0]]))
            m6_vals.append(float(row[m6_col[0]]))
    
    if len(pre_vals) < 3:
        continue
    
    mean_pre = np.mean(pre_vals)
    mean_m3 = np.mean(m3_vals)
    mean_m6 = np.mean(m6_vals)
    
    fc_m3 = mean_m3 / mean_pre if mean_pre > 0 else np.nan
    fc_m6 = mean_m6 / mean_pre if mean_pre > 0 else np.nan
    
    print(f"{g:<15} PRE={mean_pre:<8.2f} 3M={mean_m3:<8.2f}(FC={fc_m3:.3f}) 6M={mean_m6:<8.2f}(FC={fc_m6:.3f})")

# Save as JSON
with open(r'C:\Users\javie\.openclaw\workspace\lumus\belimumab_fc.json', 'w') as fp:
    json.dump(fc_data, fp, indent=2)

# Generate the report
print("\n\nGenerating report...")
lines = []
lines.append("# Belimumab RNA-seq Validation: 44 SLE Patients Pre/Post Treatment")
lines.append("")
lines.append(f"**Dataset:** Zenodo Record 14557188 (Frontiers in Immunology, 2025)")
lines.append(f"**DOI:** 10.5281/zenodo.14557188")
lines.append(f"**PMID:** 39975549")
lines.append("")
lines.append("## Study Design")
lines.append("- **SLE patients:** 44 (before treatment)")
lines.append("- **Healthy controls:** 17")
lines.append("- **Timepoints:** Pre-treatment, ~3 months post, ~6 months post")
lines.append("- **Platform:** Total RNA-seq (NovaSeq6000, PAXgene tubes)")
lines.append("- **Data:** Normalized counts (STAR alignment to hg38)")
lines.append(f"- Patients with all 3 timepoints: {len(common_ids)}")
lines.append("")
lines.append("## Results: Fold Changes (Group Means)")
lines.append("")
lines.append("| Gene | SLE vs HC | 3M vs PRE | 6M vs PRE | 6M vs HC | Interpretation |")
lines.append("|------|-----------|-----------|-----------|----------|----------------|")
for fc in fc_data:
    g = fc['Gene']
    sle = fc['SLE_vs_HC']
    m3 = fc['3M_vs_PRE']
    m6 = fc['6M_vs_PRE']
    m6hc = fc['6M_vs_HC']
    d = fc['Direction']
    lines.append(f"| {g} | {sle} | {m3} | {m6} | {m6hc} | {d} |")

lines.append("")
lines.append("## Detailed Means")
lines.append("")
lines.append("| Gene | HC | PRE | 3M POST | 6M POST |")
lines.append("|------|----|-----|---------|---------|")
for fc in fc_data:
    lines.append(f"| {fc['Gene']} | {fc['Mean_HC']} | {fc['Mean_PRE_BELIMUMAB']} | {fc['Mean_3M_POST']} | {fc['Mean_6M_POST']} |")

lines.append("")
lines.append("## Paired Analysis (Same Patients, All 3 Timepoints)")
lines.append("")
lines.append("| Gene | PRE | 3M (FC) | 6M (FC) | Direction |")
lines.append("|------|-----|---------|---------|-----------|")
for g in target_genes:
    matches = df[df['Symbol'] == g]
    if len(matches) == 0:
        continue
    row = matches.iloc[0]
    pre_vals = []
    m3_vals = []
    m6_vals = []
    for pid in common_ids:
        pre_c = [c for c in columns if f'_{pid}' in c and '_0M_' in c]
        m3_c = [c for c in columns if f'_{pid}' in c and '_3M_' in c]
        m6_c = [c for c in columns if f'_{pid}' in c and '_6M_' in c]
        if pre_c and m3_c and m6_c:
            pre_vals.append(float(row[pre_c[0]]))
            m3_vals.append(float(row[m3_c[0]]))
            m6_vals.append(float(row[m6_c[0]]))
    if len(pre_vals) < 3:
        continue
    mp = np.mean(pre_vals)
    mm3 = np.mean(m3_vals)
    mm6 = np.mean(m6_vals)
    f3 = mm3/mp if mp>0 else 0
    f6 = mm6/mp if mp>0 else 0
    direction = "BELI UP" if f6 > 1.15 else ("BELI DOWN" if f6 < 0.85 else "flat")
    lines.append(f"| {g} | {mp:.2f} | {mm3:.2f} ({f3:.3f}) | {mm6:.2f} ({f6:.3f}) | {direction} |")

lines.append("")
lines.append("## Model Validation")
lines.append("")
lines.append("**Prediction:** TLR7↑, MYD88↑, BAFF↑, IRF7↑ in SLE → should DOWNREGULATE after Belimumab")
lines.append("")
# Check key loop genes
loop_genes = ['TLR7', 'MYD88', 'TNFSF13B', 'IRF7']
for g in loop_genes:
    fc_entry = [f for f in fc_data if f['Gene'] == g]
    if fc_entry:
        f = fc_entry[0]
        sle_fc = f['SLE_vs_HC']
        beli_fc = f['6M_vs_PRE']
        sle_up = sle_fc > 1.15
        beli_down = beli_fc < 0.90
        confirmed = sle_up and beli_down
        lines.append(f"- **{g}**: SLE↑{sle_fc}x → 6M post-Belimumab: {beli_fc}x {'✅ CONFIRMED' if confirmed else '⚠️ Mixed'}")

lines.append("")
lines.append("## Limitations")
lines.append("- Normalized counts, not raw counts (can't run DESeq2)")
lines.append("- No responder/non-responder stratification available from column names")
lines.append("- Some patients missing at specific timepoints")
lines.append("- Paired analysis limited to patients with all 3 visits")
lines.append("- Multiple testing not applied at gene level")

with open(out, 'w', encoding='utf-8') as fp:
    fp.write('\n'.join(lines))

print(f"\nReport saved to {out}")
