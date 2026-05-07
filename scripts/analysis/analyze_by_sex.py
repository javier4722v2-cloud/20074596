# -*- coding: utf-8 -*-
"""
LUMUS — ANÁLISIS GTEx POR SEXO (sample-level)
TLR7, ESR1, ESR2, CYP19A1, MYD88, IRF5, TNF, TNFSF13B
"""
import sys, json, urllib.request, time
sys.stdout.reconfigure(encoding='utf-8')

API = "https://gtexportal.org/api/v2"
DS = "gtex_v8"

GENES = {
    "TLR7": "ENSG00000196664.4",
    "ESR1": "ENSG00000091831.22",
    "ESR2": "ENSG00000140009.18",
    "CYP19A1": "ENSG00000137869.14",
    "MYD88": "ENSG00000172936.12",
    "IRF5": "ENSG00000128604.19",
    "TNF": "ENSG00000232810.3",
    "TNFSF13B": "ENSG00000102524.11",
}

TISSUES = ["Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes",
           "Lung", "Skin_Sun_Exposed_Lower_leg",
           "Breast_Mammary_Tissue", "Ovary", "Uterus", "Vagina",
           "Adipose_Subcutaneous", "Liver"]

def fetch(url_path):
    url = API + url_path
    for attempt in range(3):
        try:
            req = urllib.request.Request(url)
            req.add_header("Accept", "application/json")
            resp = urllib.request.urlopen(req, timeout=30)
            return json.loads(resp.read())
        except Exception as e:
            if attempt < 2:
                time.sleep(1)
            else:
                raise

# Get ALL samples with sex for each tissue using pagination
print("🔬 LUMUS — ANÁLISIS SEXUAL GTEx v8 (sample-level)\n")

results = {}

for tid in TISSUES:
    tname = tid.replace("_", " ")
    print(f"📊 {tname}...")
    
    # Paginate through all samples
    all_samples = []
    page = 0
    while True:
        data = fetch(f"/dataset/sample?datasetId={DS}&tissueSiteDetailId={tid}&page={page}&itemsPerPage=250")
        if isinstance(data, dict) and "data" in data:
            chunk = data["data"]
            if not chunk:
                break
            all_samples.extend(chunk)
            page += 1
        else:
            break
    
    if not all_samples:
        print(f"  ⚠️ Sin datos de muestras")
        continue
    
    # Extract sample sex
    sample_sex = {}
    for s in all_samples:
        sid = s.get("sampleId", "")
        sex = s.get("sex", "")
        if sex:
            sample_sex[sid] = sex
    
    n_male = sum(1 for v in sample_sex.values() if v == "male")
    n_female = sum(1 for v in sample_sex.values() if v == "female")
    print(f"  📍 {len(sample_sex)} muestras ({n_male} H, {n_female} M)")
    
    if n_male == 0 or n_female == 0:
        print(f"  ⚠️ Solo un sexo presente, saltando")
        continue
    
    # Get expression for all genes in this tissue
    gids = ";".join(GENES.values())
    
    for sym, gid in GENES.items():
        # Initialize tissue dict for this gene if needed
        if sym not in results:
            results[sym] = {}
        if tname not in results[sym]:
            results[sym][tname] = {"M": [], "F": [], "M_n": 0, "F_n": 0}
        
        try:
            expr_data = fetch(f"/expression/geneExpression?gencodeId={gid}&datasetId={DS}&tissueSiteDetailId={tid}")
            
            if isinstance(expr_data, dict) and expr_data.get("data"):
                # The expression endpoint returns a list of dicts with "data" being the values array
                for entry in expr_data["data"]:
                    if isinstance(entry, dict) and "data" in entry:
                        values = entry["data"]
                    elif isinstance(entry, list):
                        values = entry
                    else:
                        continue
                    
                    # Map values to samples by position
                    sample_ids = list(sample_sex.keys())
                    for i, val in enumerate(values):
                        if i >= len(sample_ids):
                            break
                        sid = sample_ids[i]
                        sex = sample_sex.get(sid)
                        if sex == "male":
                            results[sym][tname]["M"].append(float(val))
                        elif sex == "female":
                            results[sym][tname]["F"].append(float(val))
                    
                    break  # Only first entry
                    
        except Exception as e:
            print(f"  ⚠️ {sym}: error getting expression: {str(e)[:50]}")
    
    print(f"  ✅ {len(GENES)} genes")

# COMPUTE and display
print(f"\n{'='*65}")
print("RESULTADOS: EXPRESIÓN PROMEDIO POR SEXO (TPM)")
print(f"{'='*65}\n")

for sym in GENES:
    print(f"\n  **{sym}**")
    for tname in [t.replace("_", " ") for t in TISSUES]:
        if tname not in results.get(sym, {}):
            continue
        td = results[sym][tname]
        m_vals = td["M"]
        f_vals = td["F"]
        
        if not m_vals or not f_vals:
            continue
        
        m_mean = sum(m_vals) / len(m_vals)
        f_mean = sum(f_vals) / len(f_vals)
        m_n = len(m_vals)
        f_n = len(f_vals)
        
        # compute std
        m_std = (sum((v-m_mean)**2 for v in m_vals) / (len(m_vals)-1))**0.5 if len(m_vals) > 1 else 0
        f_std = (sum((v-f_mean)**2 for v in f_vals) / (len(f_vals)-1))**0.5 if len(f_vals) > 1 else 0
        
        fc = f_mean / m_mean if m_mean > 0 else float('inf')
        
        results[sym][tname]["M_mean"] = round(m_mean, 3)
        results[sym][tname]["F_mean"] = round(f_mean, 3)
        results[sym][tname]["M_std"] = round(m_std, 3)
        results[sym][tname]["F_std"] = round(f_std, 3)
        results[sym][tname]["FC"] = round(fc, 3)
        results[sym][tname]["M_n"] = m_n
        results[sym][tname]["F_n"] = f_n
        
        arrow = "⬆" if fc > 1.05 else "⬇" if fc < 0.95 else "="
        print(f"    {tname:30s} M={m_mean:.3f}±{m_std:.2f}(n={m_n})  F={f_mean:.3f}±{f_std:.2f}(n={f_n})  {arrow}FC={fc:.3f}")

# Save
clean_results = {}
for sym in GENES:
    clean_results[sym] = {}
    for tname in results.get(sym, {}):
        clean_results[sym][tname] = {k: v for k, v in results[sym][tname].items() 
                                     if k not in ("M", "F")}  # Don't save raw arrays

outpath = "C:/Users/javie/.openclaw/workspace/lumus/gtex_sample_by_sex.json"
with open(outpath, "w", encoding="utf-8") as f:
    json.dump(clean_results, f, indent=2, ensure_ascii=False)
print(f"\n📁 Guardado en {outpath}")
print("✅ Completado.")
