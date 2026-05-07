import gzip, csv, json, math, re
import numpy as np
from scipy import stats as sp_stats

SOFT_FILE = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\GSE11909_family.soft.gz'
OUTPUT_CSV = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\full_de.csv'
OUTPUT_SIG = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\significant_genes.json'
ANNOT_OUT = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse11909\probe_gene_map.json'

# Step 1: Parse the SOFT file to get expression matrix and sample annotations
print("Parsing SOFT file...")
with gzip.open(SOFT_FILE, 'rt', encoding='latin-1') as f:
    content = f.read()

# Split into sections by ^
sections = content.split('^')
print(f"Found {len(sections)} sections")

# Find the series section and sample sections
series_section = None
sample_sections = []
platform_section = None

for sec in sections:
    if sec.startswith('SERIES ='):
        series_section = sec
    elif sec.startswith('SAMPLE ='):
        sample_sections.append(sec)
    elif sec.startswith('PLATFORM ='):
        platform_section = sec

print(f"Series: {series_section is not None}")
print(f"Samples: {len(sample_sections)}")
print(f"Platform: {platform_section is not None}")

# Parse sample metadata from series section
sample_metadata = {}
if series_section:
    for line in series_section.split('\n'):
        if line.startswith('!Sample_'):
            parts = line.split('=', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                val = parts[1].strip()
                if key not in sample_metadata:
                    sample_metadata[key] = []
                sample_metadata[key].append(val)

# Print sample info
for key in ['!Sample_title', '!Sample_geo_accession', '!Sample_characteristics_disease']:
    if key in sample_metadata:
        vals = sample_metadata[key]
        print(f"{key}: {len(vals)} values")
        for v in vals[:5]:
            print(f"  {v}")

# Parse individual samples from sample sections
samples = []
for sec in sample_sections:
    sample = {}
    lines = sec.split('\n')
    # Parse metadata
    i = 0
    while i < len(lines) and lines[i].startswith('!'):
        if '=' in lines[i]:
            k, v = lines[i].split('=', 1)
            sample[k.strip()] = v.strip()
        i += 1
    
    # Find table start
    while i < len(lines) and 'table_begin' not in lines[i]:
        i += 1
    
    if i < len(lines) and 'table_begin' in lines[i]:
        i += 1  # skip table_begin
        # Parse expression data
        id_vals = {}
        while i < len(lines) and 'table_end' not in lines[i]:
            parts = lines[i].strip().split('\t')
            if len(parts) >= 2:
                try:
                    id_vals[parts[0]] = float(parts[1])
                except:
                    pass
            i += 1
        sample['expression'] = id_vals
    
    samples.append(sample)

print(f"\nParsed {len(samples)} samples with expression data")

# Alternative: Try getting sample info from the SOFT file directly
# Let me try a different parsing approach - extract sample info from header
sample_info = {}
for line in series_section.split('\n') if series_section else []:
    if line.startswith('!Sample_title') or line.startswith('!Sample_description') or \
       line.startswith('!Sample_characteristics') or line.startswith('!Sample_source'):
        parts = line.split('=', 1)
        if len(parts) == 2:
            key = parts[0].strip()
            vals = [v.strip().strip('"') for v in parts[1].split('\t')]
            if key not in sample_info:
                sample_info[key] = vals
            else:
                # Extend if already has values
                sample_info[key].extend(vals)

for key, vals in sample_info.items():
    print(f"\n{key}: {len(vals)} values")
    for v in vals[:3]:
        print(f"  {v}")
