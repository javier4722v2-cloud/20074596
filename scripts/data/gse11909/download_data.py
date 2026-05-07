import GEOparse
import os

workdir = r"C:\Users\javie\.openclaw\workspace\lumus\data\gse11909"

# Try to get GSE11909 from GEO
gse = GEOparse.get_GEO(geo="GSE11909", destdir=workdir, silent=True)

# Check the platforms
print("Platforms:", list(gse.gpls.keys()))
print("Samples:", len(gse.gsms))

# Show first few sample info
for gsm_id in list(gse.gsms.keys())[:5]:
    gsm = gse.gsms[gsm_id]
    title = gsm.metadata.get('title', [''])[0]
    source = gsm.metadata.get('source_name_ch1', [''])[0]
    chars = gsm.metadata.get('characteristics_ch1', [''])
    print(f"{gsm_id}: {title} | source='{source}' | chars={chars}")

print("\n--- All sample characteristics ---")
for gsm_id in gse.gsms.keys():
    gsm = gse.gsms[gsm_id]
    chars = gsm.metadata.get('characteristics_ch1', [])
    title = gsm.metadata.get('title', [''])[0]
    source = gsm.metadata.get('source_name_ch1', [''])[0]
    geo_acc = gsm.metadata.get('geo_accession', [''])[0]
    print(f"{geo_acc}: title='{title}' | source='{source}' | chars={chars}")
