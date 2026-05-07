import gzip, csv, json, math, os
import numpy as np
from scipy import stats as sp_stats

SOFT_FILE = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\GSE11909_family.soft.gz'
OUTPUT_CSV = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\full_de.csv'
OUTPUT_SIG = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\significant_genes.json'
ANNOT_FILE = r'C:\Users\javie\.openclaw\workspace\lumus\data\GPL13158.annot.gz'  # Actually wrong platform

print(f"File size: {os.path.getsize(SOFT_FILE)} bytes")

# Read SOFT file line by line to extract sample info
sample_ids = []
sample_groups = {}  # sample_id -> 'SLE' or 'Healthy'
expression_data = {}  # probe_id -> {sample_id: value}

print("Scanning SOFT file for sample metadata...")
with gzip.open(SOFT_FILE, 'rt', encoding='latin-1') as f:
    current_sample = None
    in_table = False
    probe_col = 0
    value_col = 1
    
    for line in f:
        line = line.rstrip('\n\r')
        
        # Detect sample sections
        if line.startswith('^SAMPLE = '):
            current_sample = line.split('=', 1)[1].strip()
            sample_ids.append(current_sample)
            in_table = False
            continue
        
        # Parse metadata within sample
        if current_sample and line.startswith('!Sample_title'):
            title = line.split('=', 1)[1].strip().strip('"')
            if 'sle' in title.lower() or 'lupus' in title.lower():
                sample_groups[current_sample] = 'SLE'
            elif 'health' in title.lower() or 'normal' in title.lower() or 'control' in title.lower():
                sample_groups[current_sample] = 'Healthy'
        
        if current_sample and line.startswith('!Sample_characteristics'):
            char = line.split('=', 1)[1].strip().strip('"')
            if 'sle' in char.lower() or 'lupus' in char.lower():
                sample_groups[current_sample] = 'SLE'
            elif 'health' in char.lower() or 'normal' in char.lower() or 'control' in char.lower():
                sample_groups[current_sample] = 'Healthy'
        
        # Detect table start
        if 'table_begin' in line:
            in_table = True
            continue
        
        if 'table_end' in line:
            in_table = False
            current_sample = None
            continue
        
        # Parse expression
        if in_table and current_sample:
            parts = line.split('\t')
            if len(parts) >= 2:
                probe = parts[0].strip()
                try:
                    val = float(parts[1].strip())
                except:
                    continue
                
                if probe not in expression_data:
                    expression_data[probe] = {}
                expression_data[probe][current_sample] = val

print(f"Found {len(sample_ids)} samples")
sle_samples = [s for s, g in sample_groups.items() if g == 'SLE']
healthy_samples = [s for s, g in sample_groups.items() if g == 'Healthy']
print(f"SLE: {len(sle_samples)}, Healthy: {len(healthy_samples)}, Unknown: {len(sample_ids) - len(sle_samples) - len(healthy_samples)}")

if len(sle_samples) == 0:
    print("No SLE/Healthy labels found. Printing sample titles...")
    # Try reading titles differently
    with gzip.open(SOFT_FILE, 'rt', encoding='latin-1') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('!Sample_title'):
                parts = line.split('\t')
                titles = [p.strip().strip('"') for p in parts[1:]]
                for i, t in enumerate(titles):
                    sid = sample_ids[i] if i < len(sample_ids) else f'GSM{i}'
                    if 'sle' in t.lower():
                        sample_groups[sid] = 'SLE'
                    elif 'health' in t.lower() or 'normal' in t.lower():
                        sample_groups[sid] = 'Healthy'
                break
    
    sle_samples = [s for s, g in sample_groups.items() if g == 'SLE']
    healthy_samples = [s for s, g in sample_groups.items() if g == 'Healthy']
    print(f"After retry: SLE: {len(sle_samples)}, Healthy: {len(healthy_samples)}")

print(f"\nProbes collected: {len(expression_data)}")
print(f"Example SLE: {sle_samples[:3]}")
print(f"Example Healthy: {healthy_samples[:3]}")

# Perform DE
print("\nPerforming differential expression...")
results = []
for probe, vals in expression_data.items():
    sle_vals = [v for s, v in vals.items() if s in sle_samples]
    h_vals = [v for s, v in vals.items() if s in healthy_samples]
    
    if len(sle_vals) < 3 or len(h_vals) < 3:
        continue
    
    sle_mean = np.mean(sle_vals)
    h_mean = np.mean(h_vals)
    
    if sle_mean == 0 and h_mean == 0:
        continue
    
    fc = sle_mean / h_mean if h_mean > 0 else 1.0
    log2fc = math.log2(max(fc, 0.001))
    
    try:
        t_stat, p_val = sp_stats.ttest_ind(sle_vals, h_vals, equal_var=False)
    except:
        p_val = 1.0
    
    results.append({
        'probe': probe,
        'sle_mean': round(sle_mean, 4),
        'healthy_mean': round(h_mean, 4),
        'fc': round(fc, 4),
        'log2fc': round(log2fc, 4),
        'pval': p_val,
        'n_sle': len(sle_vals),
        'n_healthy': len(h_vals)
    })

print(f"Processed {len(results)} probes")

# Save
print("Saving...")
with open(OUTPUT_CSV, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['probe', 'sle_mean', 'healthy_mean', 'fc', 'log2fc', 'pval', 'n_sle', 'n_healthy'])
    w.writeheader()
    w.writerows(results)

sig = [r for r in results if r['pval'] < 0.05]
sig_up = [r for r in sig if r['log2fc'] > 0.5]
sig_down = [r for r in sig if r['log2fc'] < -0.5]

print(f"\n=== RESULTS GSE11909 ===")
print(f"Probes analyzed: {len(results)}")
print(f"Significant (p<0.05): {len(sig)}")
print(f"UP >0.5: {len(sig_up)}")
print(f"DOWN <-0.5: {len(sig_down)}")

print(f"\nTop 20 UP:")
for r in sorted(sig, key=lambda x: -x['log2fc'])[:20]:
    print(f"  {r['probe']}: FC={r['fc']:.2f}, log2FC={r['log2fc']:.2f}, p={r['pval']:.2e}")

print(f"\nTop 20 DOWN:")
for r in sorted(sig, key=lambda x: x['log2fc'])[:20]:
    print(f"  {r['probe']}: FC={r['fc']:.2f}, log2FC={r['log2fc']:.2f}, p={r['pval']:.2e}")

# Now download GPL96 annotation and map probes to genes
print("\nDownloading GPL96 annotation...")
import urllib.request
url = 'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL96nnn/GPL96/annot/GPL96.annot.gz'
annot_local = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\GPL96.annot.gz'
try:
    urllib.request.urlretrieve(url, annot_local)
    print("Downloaded GPL96 annotation")
    
    # Load probe->gene map
    probe_gene = {}
    with gzip.open(annot_local, 'rt', encoding='latin-1') as f:
        for line in f:
            if line.startswith('ID\t') or line.startswith('!'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                pid = parts[0].strip()
                gene = parts[6].strip()
                if gene and gene != '---':
                    probe_gene[pid] = gene
    
    print(f"Loaded {len(probe_gene)} probe->gene mappings")
    
    # Map results to genes
    print(f"\n=== GENE-MAPPED RESULTS ===")
    for r in sorted(sig, key=lambda x: -abs(x['log2fc']))[:30]:
        gene = probe_gene.get(r['probe'], r['probe'])
        d = "UP" if r['log2fc'] > 0 else "DOWN"
        print(f"  {gene:15s} {d}: FC={r['fc']:.2f}, log2FC={r['log2fc']:.2f}, p={r['pval']:.2e}")
    
    # Check loop genes
    loop_genes = {'TLR7','TLR8','MYD88','IRF7','IFNA1','TNFSF13B','TNFRSF13C','NFKB1','NFKB2',
                  'CXCL10','TNF','CYP19A1','ESR1','STAT1','STAT3','IL6','IL10','CD19','MS4A1'}
    print(f"\n=== LOOP GENES IN GSE11909 ===")
    found = []
    for r in results:
        gene = probe_gene.get(r['probe'], '')
        if gene in loop_genes:
            d = "UP" if r['log2fc'] > 0 else "DOWN"
            sig_mark = "*" if r['pval'] < 0.05 else ""
            found.append(gene)
            print(f"  {gene:15s} {d}{sig_mark}: FC={r['fc']:.2f}, log2FC={r['log2fc']:.2f}, p={r['pval']:.2e}")
    missing = loop_genes - set(found)
    print(f"Not on GPL96 chip: {missing}")
    
except Exception as e:
    print(f"Could not download annotation: {e}")
    print("Showing probe-level results (no gene mapping)")
