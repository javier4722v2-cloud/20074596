"""Analisis del loop TLR7/BAFF/IFN en SLE adulto (GSE81622)."""
import numpy as np
from scipy.stats import pearsonr, spearmanr
import csv

def load_csv(path):
    with open(path) as f:
        reader = csv.DictReader(f)
        data = list(reader)
    genes = list(data[0].keys())
    matrix = {g: [] for g in genes}
    for row in data:
        for g in genes:
            matrix[g].append(float(row[g]))
    return genes, {g: np.array(v) for g, v in matrix.items()}

# Cargar datos
sle_genes, sle = load_csv(r'C:\Users\javie\.openclaw\workspace\lumus\loop_analysis\SLE_adult\SLE.csv')
healthy_genes, healthy = load_csv(r'C:\Users\javie\.openclaw\workspace\lumus\loop_analysis\SLE_adult\Healthy.csv')

assert sle_genes == healthy_genes, "Gene mismatch!"
genes = sle_genes

# --- 1. Medias y fold changes ---
results = []
for g in genes:
    m_sle = np.mean(sle[g])
    m_h = np.mean(healthy[g])
    fc = m_sle / m_h if m_h != 0 else float('inf')
    log2fc = np.log2(fc) if fc > 0 else float('inf')
    results.append((g, m_sle, m_h, fc, log2fc))

results.sort(key=lambda x: x[3], reverse=True)

# --- 2. Correlaciones entre pares clave ---
pairs = [
    ("TLR7", "IRF7"),
    ("IRF7", "IFIT1"),
    ("TNFSF13B", "IRF7"),
    ("IFIT1", "ISG15"),
    ("IFIT1", "RSAD2"),
]

corr_results = []
for g1, g2 in pairs:
    if g1 not in sle or g2 not in sle:
        continue
    r_sle, p_sle = pearsonr(sle[g1], sle[g2])
    rho_sle, p_rho_sle = spearmanr(sle[g1], sle[g2])
    r_h, p_h = pearsonr(healthy[g1], healthy[g2])
    rho_h, p_rho_h = spearmanr(healthy[g1], healthy[g2])
    corr_results.append((g1, g2, r_sle, p_sle, rho_sle, r_h, p_h, rho_h))

# --- Construir output ---
lines = []
lines.append("=" * 80)
lines.append("ANALISIS DEL LOOP TLR7/BAFF/IFN - SLE ADULTO (GSE81622)")
lines.append("30 SLE adultas vs 25 Healthy adultas")
lines.append("=" * 80)

lines.append("\n\n1. GENES DEL LOOP - MEDIAS Y FOLD CHANGES")
lines.append("-" * 80)
lines.append(f"{'Gen':<12} {'SLE Media':<14} {'Healthy Media':<14} {'Fold Chg':<14} {'log2FC':<10}")
lines.append("-" * 80)

loop_genes = ["TLR7", "MYD88", "IRAK4", "TICAM1", "TNFSF13B", "CD40LG",
               "IRF7", "STAT1", "STAT3", "IFIT1", "IFIT3", "ISG15", "IFI44",
               "IFI44L", "IFI6", "MX1", "OAS1", "RSAD2", "HERC5", "DDX58",
               "IFIH1", "CXCL10", "CCL5", "TNF", "IL1B", "IL6", "IFNA1", "IFI27"]

for g, m_sle, m_h, fc, l2fc in results:
    if g in loop_genes:
        lines.append(f"{g:<12} {m_sle:<14.2f} {m_h:<14.2f} {fc:<14.4f} {l2fc:<+10.4f}")

lines.append("\n\nTODOS LOS GENES ORDENADOS POR FOLD CHANGE")
lines.append("-" * 80)
lines.append(f"{'Gen':<12} {'SLE Media':<14} {'Healthy Media':<14} {'Fold Chg':<14} {'log2FC':<10}")
lines.append("-" * 80)
for g, m_sle, m_h, fc, l2fc in results:
    lines.append(f"{g:<12} {m_sle:<14.2f} {m_h:<14.2f} {fc:<14.4f} {l2fc:<+10.4f}")

lines.append("\n\n2. CORRELACIONES (PEARSON) - PARES CLAVE DEL LOOP")
lines.append("-" * 80)
lines.append(f"{'Par':<25} {'SLE r':<10} {'SLE p-val':<12} {'SLE rho':<10} {'Healthy r':<10} {'Healthy p-val':<12} {'Healthy rho':<8}")
lines.append("-" * 80)

for g1, g2, r_sle, p_sle, rho_sle, r_h, p_h, rho_h in corr_results:
    label = f"{g1}->{g2}"
    lines.append(f"{label:<25} {r_sle:<+10.4f} {p_sle:<12.2e} {rho_sle:<+10.4f} {r_h:<+10.4f} {p_h:<12.2e} {rho_h:<+8.4f}")

lines.append("\n\n3. EVALUACION DEL LOOP (comparacion con SLE pediatrico)")
lines.append("-" * 80)

# Analisis de activacion del loop
lines.append("\n--- SE~NALES DE ACTIVACION DEL LOOP ---")
fc_map = {r[0]: r[3] for r in results}

key_genes_status = [
    ("TLR7", fc_map["TLR7"]),
    ("IRF7", fc_map["IRF7"]),
    ("TNFSF13B (BAFF)", fc_map["TNFSF13B"]),
    ("IFIT1", fc_map["IFIT1"]),
    ("ISG15", fc_map["ISG15"]),
    ("RSAD2", fc_map["RSAD2"]),
    ("IFNA1", fc_map["IFNA1"]),
    ("IFI27", fc_map["IFI27"]),
]

for g, fc_val in key_genes_status:
    if fc_val > 1:
        lines.append(f"  [+] {g}: UP {fc_val:.2f}x sobreexpresado en SLE adulto")
    elif fc_val < 1:
        lines.append(f"  [-] {g}: DOWN {fc_val:.2f}x reprimido en SLE adulto")
    else:
        lines.append(f"  [~] {g}: {fc_val:.2f}x sin cambio significativo")

# Fold change del IFNA1 -> IFI27 (IFN signature)
fc_ifna1 = fc_map["IFNA1"]
fc_ifi27 = fc_map["IFI27"]
fc_mx1 = fc_map["MX1"]
fc_cxcl10 = fc_map["CXCL10"]
fc_ccr5 = fc_map.get("CCL5", 0)

lines.append(f"\n  Firma IFN tipo I: IFNA1={fc_ifna1:.2f}x, IFI27={fc_ifi27:.2f}x, MX1={fc_mx1:.2f}x, CXCL10={fc_cxcl10:.2f}x")

lines.append("\n--- CORRELACIONES vs UMBRAL PEDIATRICO (r > 0.7 esperado) ---")

pediatric_ref = {
    ("TLR7", "IRF7"): "> 0.7 (activacion del eje TLR7-IRF7)",
    ("IRF7", "IFIT1"): "> 0.7 (respuesta IFN corriente abajo)",
    ("TNFSF13B", "IRF7"): "> 0.7 (BAFF activando via IRF7)",
    ("IFIT1", "ISG15"): "> 0.7 (genes estimulados por IFN co-regulados)",
    ("IFIT1", "RSAD2"): "> 0.7 (firma IFN consolidada)",
}

for g1, g2, r_sle, p_sle, rho_sle, r_h, p_h, rho_h in corr_results:
    label = f"{g1}<->{g2}"
    expected = pediatric_ref.get((g1, g2), "desconocido")

    if r_sle > 0.7:
        status = "[OK] CONFIRMADO"
        note = f"(r={r_sle:.3f} > 0.7, similar a pediatrico)"
    elif r_sle > 0.5:
        status = "[~] PARCIAL"
        note = f"(r={r_sle:.3f}, moderada pero < 0.7)"
    else:
        status = "[X] NO"
        note = f"(r={r_sle:.3f}, no alcanza umbral)"

    lines.append(f"  {status}: {label} -> r={r_sle:.4f} {note}")
    lines.append(f"         Esperado pediatrico: {expected}")

    if abs(r_sle - r_h) > 0.3:
        lines.append(f"         Diferencia marcada con Healthy (r_Healthy={r_h:.4f})")

lines.append("\n\n--- RESUMEN ---")
strong_corrs = sum(1 for _, _, r, _, _, _, _, _ in corr_results if r > 0.7)
mod_corrs = sum(1 for _, _, r, _, _, _, _, _ in corr_results if r > 0.5 and r <= 0.7)
weak_corrs = sum(1 for _, _, r, _, _, _, _, _ in corr_results if r <= 0.5)

lines.append(f"  Correlaciones fuertes (> 0.7): {strong_corrs}/{len(corr_results)}")
lines.append(f"  Correlaciones moderadas (0.5-0.7): {mod_corrs}/{len(corr_results)}")
lines.append(f"  Correlaciones debiles (<= 0.5): {weak_corrs}/{len(corr_results)}")

# Genes IFN estimulados upregulados
ifn_genes = ["IFIT1","IFIT3","ISG15","MX1","OAS1","OAS2","RSAD2","IFI44","IFI44L","IFI6","HERC5"]
ifn_up = sum(1 for g, _, _, fc, _ in results if g in ifn_genes and fc > 1.5)
lines.append(f"  Genes IFN-estimulados con FC > 1.5: {ifn_up}/{len(ifn_genes)}")

fc_avg_loop = np.mean([fc_map["TLR7"], fc_map["IRF7"], fc_map["TNFSF13B"],
                        fc_map["IFIT1"], fc_map["ISG15"], fc_map["RSAD2"]])
lines.append(f"  Fold change promedio del loop: {fc_avg_loop:.3f}x")

# Correlaciones Spearman como control de robustez
lines.append("\n--- CORRELACIONES SPEARMAN (control de robustez) ---")
lines.append(f"{'Par':<25} {'SLE rho':<10} {'Healthy rho':<10}")
lines.append("-" * 50)
for g1, g2, r_sle, p_sle, rho_sle, r_h, p_h, rho_h in corr_results:
    label = f"{g1}<->{g2}"
    lines.append(f"{label:<25} {rho_sle:<+10.4f} {rho_h:<+10.4f}")

if fc_avg_loop > 1.5 and strong_corrs >= 3:
    conclusion = ("CONCLUSION: El loop TLR7/BAFF/IFN esta ACTIVO en SLE adulto "
                  "(GSE81622), con sobreexpresion generalizada y correlaciones "
                  "que confirman la senalizacion en cascada.")
elif fc_avg_loop > 1.2 and strong_corrs >= 2:
    conclusion = ("CONCLUSION: El loop TLR7/BAFF/IFN esta PARCIALMENTE activo. "
                  "Hay sobreexpresion de varios genes pero las correlaciones no "
                  "alcanzan el nivel del SLE pediatrico.")
else:
    conclusion = ("CONCLUSION: El loop TLR7/BAFF/IFN NO se confirma claramente "
                  "en SLE adulto (GSE81622). Las correlaciones y/o sobreexpresion "
                  "son debiles.")

lines.append(f"\n  {conclusion}")

lines.append("\n\n" + "=" * 80)
lines.append("FIN DEL ANALISIS")
lines.append("=" * 80)

output = "\n".join(lines)

# Guardar
with open(r'C:\Users\javie\.openclaw\workspace\lumus\loop_analysis\SLE_adult\analysis.txt', 'w', encoding='utf-8') as f:
    f.write(output)

print("[OK] analysis.txt guardado")
print("=" * 80)
# Mostrar resumen rapido
for line in output.split("\n"):
    if any(kw in line for kw in ["GENES DEL LOOP", "CORRELACIONES (PEARSON)", "RESUMEN", "CONCLUSION", "CONFIRMADO", "PARCIAL", "NO ", "ACTIVO", "promedio"]):
        print(line)
