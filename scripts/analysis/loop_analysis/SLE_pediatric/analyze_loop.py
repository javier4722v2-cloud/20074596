import pandas as pd
import numpy as np

sle = pd.read_csv('SLE.csv')
healthy = pd.read_csv('Healthy.csv')

genes = ['TLR7','MYD88','IRAK4','TICAM1','TNFSF13B','CD40LG','IRF7','STAT1','STAT3',
         'IFIT1','IFIT3','ISG15','IFI44','IFI44L','IFI6','MX1','OAS1','OAS2','OAS3',
         'OASL','RSAD2','HERC5','HERC6','XAF1','BATF2','EPSTI1','DDX58','DDX60',
         'IFIH1','CMPK2','IFITM1','CXCL10','CCL5','TNF','IL1B','IL6','ELANE',
         'ESR1','CYP19A1','IFNA1','IFI27']

sle_mean = sle[genes].mean()
healthy_mean = healthy[genes].mean()
fc = sle_mean - healthy_mean
fc_lin = 2 ** fc

tbl = pd.DataFrame({'SLE_mean': sle_mean, 'Healthy_mean': healthy_mean, 'FC_log2': fc, 'FC_linear': fc_lin})
tbl = tbl.sort_values('FC_linear', ascending=False)

print('=' * 80)
print('1. MEDIAS POR GEN | SLE vs HEALTHY')
print('=' * 80)
print(tbl.round(3).to_string())

pairs = [
    ('TLR7','MYD88'), ('TLR7','IRF7'), ('IRF7','IFIT1'), ('IRF7','ISG15'),
    ('TNFSF13B','IRF7'), ('IFIT1','IFI44L'), ('IFIT1','RSAD2'), ('MX1','OAS1'),
    ('RSAD2','CXCL10')
]

print()
print('=' * 80)
print('2. CORRELACIONES PARES CLAVE')
print('=' * 80)

results = []
for g1, g2 in pairs:
    rs = sle[g1].corr(sle[g2])
    rh = healthy[g1].corr(healthy[g2])
    results.append((g1, g2, rs, rh))
    ss = 'OK' if abs(rs) > 0.7 else 'DEBIL'
    sh = 'OK' if abs(rh) > 0.7 else 'DEBIL'
    line = "{:10s} -> {:8s}   SLE: r={:.4f} [{}]   Healthy: r={:.4f} [{}]".format(g1, g2, rs, ss, rh, sh)
    print(line)

print()
print('=' * 80)
print('3. CONCLUSION')
print('=' * 80)

weak = [(g1, g2, rs, rh) for g1, g2, rs, rh in results if abs(rs) < 0.7 or abs(rh) < 0.7]
if not weak:
    print('LOOP CONFIRMADO: todos los pares tienen r > 0.7 en SLE y Healthy.')
else:
    print('Pares con correlacion debil (r < 0.7):')
    for g1, g2, rs, rh in weak:
        parts = []
        if abs(rs) < 0.7:
            parts.append("SLE r={:.4f}".format(rs))
        if abs(rh) < 0.7:
            parts.append("Healthy r={:.4f}".format(rh))
        print("  {:10s} -> {:8s}: {}".format(g1, g2, ", ".join(parts)))

with open('analysis.txt', 'w', encoding='utf-8') as f:
    f.write('=== LOOP TLR7/BAFF/IFN | SLE PEDIATRICO (GSE65391) ===\n')
    f.write('Muestras: SLE={}, Healthy={}\n\n'.format(len(sle), len(healthy)))
    f.write('1. MEDIAS POR GEN\n')
    f.write(tbl.round(3).to_string())
    f.write('\n\n')
    f.write('2. CORRELACIONES PARES CLAVE\n')
    for g1, g2, rs, rh in results:
        f.write('{} -> {}:  SLE r={:.4f}  Healthy r={:.4f}\n'.format(g1, g2, rs, rh))
    f.write('\n3. CONCLUSION\n')
    if not weak:
        f.write('Loop CONFIRMADO: todas las correlaciones r > 0.7 en SLE y Healthy.\n')
    else:
        f.write('Loop NO completamente confirmado. Pares debiles:\n')
        for g1, g2, rs, rh in weak:
            parts = []
            if abs(rs) < 0.7:
                parts.append("SLE r={:.4f}".format(rs))
            if abs(rh) < 0.7:
                parts.append("Healthy r={:.4f}".format(rh))
            f.write('  {} -> {}: {}\n'.format(g1, g2, ", ".join(parts)))

print()
print('Resultados guardados en analysis.txt')
