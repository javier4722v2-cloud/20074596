"""Analyze loop by organ involvement (SLEDAI components) in GSE65391 blood"""
import gzip, csv, os, json
import numpy as np
from collections import defaultdict

DATA_DIR = r'C:\Users\javie\.openclaw\workspace\lumus'
OUT_DIR = DATA_DIR + r'\loop_analysis\by_organ'
os.makedirs(OUT_DIR, exist_ok=True)

LOOP_GENES = [
    'TLR7', 'MYD88', 'IRAK4', 'TICAM1', 'TNFSF13B', 'CD40LG',
    'IRF7', 'STAT1', 'STAT3',
    'IFIT1', 'IFIT3', 'ISG15', 'IFI44', 'IFI44L', 'IFI6', 'MX1',
    'OAS1', 'OAS2', 'OAS3', 'OASL',
    'RSAD2', 'HERC5', 'HERC6', 'XAF1', 'BATF2', 'EPSTI1',
    'DDX58', 'DDX60', 'IFIH1', 'CMPK2', 'IFITM1',
    'CXCL10', 'CCL5', 'TNF', 'IL1B', 'IL6',
    'ELANE', 'ESR1', 'CYP19A1', 'IFNA1', 'IFI27'
]

# ===== Step 1: Parse metadata =====
f = DATA_DIR + r'\cle_data\GSE65391_series_matrix.txt.gz'
lines = gzip.open(f, 'rt', encoding='utf-8', errors='replace').readlines()

titles = [s.strip().strip('"') for s in [l for l in lines if l.startswith('!Sample_title')][0].split('\t')[1:]]

# Find all characteristics lines
char_fields = {}
for l in lines:
    if l.startswith('!Sample_characteristics_ch1'):
        parts = l.split('\t')
        # First sample tells us the field name
        first_val = parts[1].strip().strip('"')
        field_name = first_val.split(':')[0].strip() if ':' in first_val else first_val
        char_fields[field_name] = [None] * len(titles)
        
        for i, p in enumerate(parts[1:]):
            if i >= len(titles): break
            val = p.strip().strip('"')
            if ':' in val:
                char_fields[field_name][i] = val.split(':', 1)[1].strip()
            else:
                char_fields[field_name][i] = val

# ===== Step 2: Define organ systems from SLEDAI components =====
# SLEDAI component scoring: 0 = absent, non-zero = present
# The number indicates SLEDAI weight, not severity
ORGAN_SYSTEMS = {
    'CNS': ['seizure', 'psychosis', 'organic_brain_syndrome', 'visual_disturbance',
            'cranial_nerve_disorder', 'lupus_headache', 'cva'],
    'Renal': ['urinary_casts', 'hematuria', 'proteinuria', 'pyuria', 'renal'],
    'Mucocutaneous': ['new_rash', 'alopecia', 'mucosal_ulcers'],
    'Musculoskeletal': ['arthritis', 'myositis'],
    'Serosal': ['pleurisy', 'pericarditis'],
    'Hematologic': ['thrombocytopenia', 'leukopenia'],
    'Vascular': ['vasculitis'],
    'Fever': ['fever'],
}

# Categorize each SLE sample by which organ systems are affected
n = len(titles)
disease = char_fields.get('disease state', [''] * n)

print(f'Total samples: {n}')
sle_samples = sum(1 for d in disease if d == 'SLE')
healthy_samples = sum(1 for d in disease if d == 'Healthy')
print(f'SLE: {sle_samples}, Healthy: {healthy_samples}')

# For each SLE sample, determine which organ systems are involved
sample_organs = {i: [] for i in range(n) if disease[i] == 'SLE'}

for org_name, components in ORGAN_SYSTEMS.items():
    for comp in components:
        values = char_fields.get(comp, [])
        if not values or len(values) != n:
            continue
        for i in list(sample_organs.keys()):
            val = values[i]
            try:
                v = float(val)
                if v > 0:  # Non-zero = organ involved
                    if org_name not in sample_organs[i]:
                        sample_organs[i].append(org_name)
            except:
                pass

# Count per organ system
organ_counts = defaultdict(int)
for organs in sample_organs.values():
    for o in organs:
        organ_counts[o] += 1

print('\nOrgan system involvement in SLE:')
for o, c in sorted(organ_counts.items(), key=lambda x: -x[1]):
    print(f'  {o:<18s} {c:>4d} samples')

# ===== Step 3: Extract expression data =====
print('\nExtracting expression data...')
probe_to_gene = {}
with open(DATA_DIR + r'\cle_data\probe_to_gene.txt', encoding='utf-8') as fh:
    for line in fh:
        p = line.strip().split('\t')
        if len(p) >= 2: probe_to_gene[p[0]] = p[1]

target_probes = set(p for p, g in probe_to_gene.items() if g in LOOP_GENES)

# Find table
ts = te = 0
for i, l in enumerate(lines):
    if 'series_matrix_table_begin' in l:
        ts = i + 1
        break
for i in range(ts, len(lines)):
    if 'series_matrix_table_end' in lines[i]:
        te = i
        break

gene_sample_vals = defaultdict(lambda: defaultdict(list))

for i in range(ts, te):
    parts = lines[i].strip().split('\t')
    probe = parts[0].strip().strip('"')
    if probe not in target_probes:
        continue
    gene = probe_to_gene[probe]
    for j in range(n):
        col = j + 1
        if col >= len(parts): continue
        try:
            v = float(parts[col].strip().strip('"'))
            gene_sample_vals[gene][j].append(v)
        except:
            pass

# Aggregate per gene per sample
expr = {g: {} for g in LOOP_GENES}
for j in range(n):
    for g in LOOP_GENES:
        vals = gene_sample_vals[g].get(j, [])
        if vals:
            expr[g][j] = float(np.mean(vals))

# ===== Step 4: Compute correlations per organ system =====
KEY_PAIRS = [
    ('TLR7', 'IRF7'), ('IRF7', 'IFIT1'), ('IRF7', 'ISG15'),
    ('TNFSF13B', 'IRF7'), ('IFIT1', 'IFI44L'), ('IFIT1', 'RSAD2'),
    ('MX1', 'OAS1'), ('RSAD2', 'CXCL10')
]

# Healthy samples
healthy_idx = [i for i in range(n) if disease[i] == 'Healthy']

def correlation(xs, ys):
    """Pearson r"""
    if len(xs) < 3: return 0, 1
    xa, ya = np.array(xs), np.array(ys)
    mx, my = np.nanmean(xa), np.nanmean(ya)
    xc, yc = xa - mx, ya - my
    num = np.nansum(xc * yc)
    den = np.sqrt(np.nansum(xc**2) * np.nansum(yc**2))
    if den == 0: return 0, 1
    r = num / den
    return r, 1 if np.isnan(r) else r

results = {}
results['Healthy'] = {}
for g1, g2 in KEY_PAIRS:
    xs = [expr[g1][i] for i in healthy_idx if i in expr[g1] and i in expr[g2]]
    ys = [expr[g2][i] for i in healthy_idx if i in expr[g1] and i in expr[g2]]
    r, _ = correlation(xs, ys)
    results['Healthy'][f'{g1}->{g2}'] = {'r': float(r), 'n': len(xs)}

# All SLE
all_sle_idx = [i for i in range(n) if disease[i] == 'SLE']
results['SLE_All'] = {}
for g1, g2 in KEY_PAIRS:
    xs = [expr[g1][i] for i in all_sle_idx if i in expr[g1] and i in expr[g2]]
    ys = [expr[g2][i] for i in all_sle_idx if i in expr[g1] and i in expr[g2]]
    r, _ = correlation(xs, ys)
    results['SLE_All'][f'{g1}->{g2}'] = {'r': float(r), 'n': len(xs)}

# Per organ system
all_org_names = sorted(organ_counts.keys())
for org in all_org_names:
    org_idx = [i for i in all_sle_idx if org in sample_organs.get(i, [])]
    # Also get "NOT this organ" for comparison
    not_org_idx = [i for i in all_sle_idx if org not in sample_organs.get(i, [])]
    
    results[org] = {}
    results[f'No_{org}'] = {}
    
    for g1, g2 in KEY_PAIRS:
        # Samples WITH this organ involvement
        xs = [expr[g1][i] for i in org_idx if i in expr[g1] and i in expr[g2]]
        ys = [expr[g2][i] for i in org_idx if i in expr[g1] and i in expr[g2]]
        r, _ = correlation(xs, ys)
        results[org][f'{g1}->{g2}'] = {'r': float(r), 'n': len(xs)}
        
        # Samples WITHOUT this organ involvement
        xs = [expr[g1][i] for i in not_org_idx if i in expr[g1] and i in expr[g2]]
        ys = [expr[g2][i] for i in not_org_idx if i in expr[g1] and i in expr[g2]]
        r, _ = correlation(xs, ys)
        results[f'No_{org}'][f'{g1}->{g2}'] = {'r': float(r), 'n': len(xs)}

# ===== Step 5: Save and report =====
with open(OUT_DIR + '/results.json', 'w') as f:
    json.dump(results, f, indent=2)

# Also compute means per organ for top loop genes
TOP_GENES = ['IFI27', 'RSAD2', 'IFI44L', 'IFIT1', 'ISG15', 'IRF7', 'TNFSF13B', 'TLR7']
mean_results = {}
for org in ['Healthy'] + ['SLE_All'] + all_org_names:
    if org == 'Healthy':
        idx = healthy_idx
    elif org == 'SLE_All':
        idx = all_sle_idx
    else:
        idx = [i for i in all_sle_idx if org in sample_organs.get(i, [])]
    
    group = {}
    for g in TOP_GENES:
        vals = [expr[g][i] for i in idx if i in expr[g]]
        if vals:
            group[g] = {
                'mean': float(np.mean(vals)),
                'n': len(vals)
            }
    mean_results[org] = group
    
    # Also get NOT this organ
    if org not in ('Healthy', 'SLE_All'):
        idx_no = [i for i in all_sle_idx if org not in sample_organs.get(i, [])]
        group_no = {}
        for g in TOP_GENES:
            vals = [expr[g][i] for i in idx_no if i in expr[g]]
            if vals:
                group_no[g] = {
                    'mean': float(np.mean(vals)),
                    'n': len(vals)
                }
        mean_results[f'No_{org}'] = group_no

with open(OUT_DIR + '/means.json', 'w') as f:
    json.dump(mean_results, f, indent=2)

# ===== Print report =====
lines_out = []
lines_out.append('=' * 70)
lines_out.append('LOOP ANALISIS POR AFECTACION ORGANICA (SLEDAI components)')
lines_out.append('Dataset: GSE65391 (pediatric SLE blood, n=924+72)')
lines_out.append('=' * 70)
lines_out.append('')
lines_out.append(f'{'Organ system':<20s}  {'n':>5s}  {'% SLE':>6s}')
lines_out.append('-' * 35)
for org in all_org_names:
    cnt = organ_counts[org]
    pct = cnt / sle_samples * 100
    lines_out.append(f'{org:<20s}  {cnt:>5d}  {pct:>5.1f}%')

lines_out.append('')
lines_out.append('')
lines_out.append('CORRELACIONES POR ORGANO AFECTADO')
lines_out.append('=' * 70)

# Header
header = f'{'Grupo':<22s}'
for g1, g2 in KEY_PAIRS:
    header += f'  {g1[:4]}->{g2[:4]}'
lines_out.append(header)
lines_out.append('-' * len(header))

def format_row(name, data_dict, pairs=KEY_PAIRS):
    row = f'{name:<22s}'
    for g1, g2 in pairs:
        key = f'{g1}->{g2}'
        if key in data_dict:
            r = data_dict[key]['r']
            n = data_dict[key]['n']
            if r > 0.7: sym = '  ✅'
            elif r > 0.5: sym = '  🟠'
            else: sym = '  ❌'
            row += f'  {sym}'
        else:
            row += f'  ？'
    return row

# First: Healthy and SLE All
lines_out.append(format_row('Healthy', results['Healthy']))
lines_out.append(format_row('SLE_All', results['SLE_All']))
lines_out.append('')

for org in all_org_names:
    lines_out.append(format_row(f'With_{org}', results[org]))
    lines_out.append(format_row(f'No_{org}', results.get(f'No_{org}', {})))

lines_out.append('')
lines_out.append('')
lines_out.append('EXPRESION MEDIA DE GENES CLAVE')
lines_out.append('=' * 70)

header2 = f'{'Grupo':<20s}'
for g in TOP_GENES:
    header2 += f'  {g:<10s}'
lines_out.append(header2)
lines_out.append('-' * len(header2))

for org in ['Healthy', 'SLE_All'] + all_org_names:
    row = f'{org:<20s}'
    for g in TOP_GENES:
        if g in mean_results.get(org, {}):
            m = mean_results[org][g]['mean']
            row += f'  {m:>8.2f} '
        else:
            row += f'  {"--":>8s} '
    lines_out.append(row)

lines_out.append('')
lines_out.append('')
lines_out.append('COMPARACION: CON vs SIN afectacion (fold change del gen respecto a Healthy)')
lines_out.append('=' * 70)

# Compute FC vs Healthy for each group
healthy_means = {}
for g in TOP_GENES:
    vals = [expr[g][i] for i in healthy_idx if i in expr[g]]
    healthy_means[g] = float(np.mean(vals)) if vals else None

for org in all_org_names:
    lines_out.append(f'\n--- {org}: FC vs Healthy ---')
    idx_with = [i for i in all_sle_idx if org in sample_organs.get(i, [])]
    idx_without = [i for i in all_sle_idx if org not in sample_organs.get(i, [])]
    
    for label, idx in [('  Con ', idx_with), ('  Sin ', idx_without)]:
        row = label
        for g in TOP_GENES:
            vals = [expr[g][i] for i in idx if i in expr[g]]
            if vals and healthy_means[g]:
                fc = float(np.mean(vals)) / healthy_means[g]
                row += f'  {g}: {fc:.2f}x'
        lines_out.append(row)

with open(OUT_DIR + '/analysis.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines_out))

print('\n'.join(lines_out))
print(f'\nResults saved to {OUT_DIR}/')
