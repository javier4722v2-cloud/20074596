"""
Healthy Expression Reference
Auto-generated from GSE65391 (pediatric) + GSE81622 (adult)
Platform: GPL10558 (Illumina HumanHT-12)

Usage:
    from healthy_ref import ref, compare_to_healthy, multi_gene_report

    ref["Pediatric"]["genes"]["TLR7"]["mean"]
    ref["Adult"]["genes"]["IRF7"]["mean"]
    ref["All"]["genes"]["BAFF"]["mean"]
"""

import numpy as np


reference = {
    "Pediatric": {
        "n": 72,
        "source": "GSE65391 (pediatric whole blood, 72 samples)",
        "genes": {
            "TLR7": {
                "mean": 6.9787,
                "std": 0.1025,
                "n": 72
            },
            "TLR8": {
                "mean": 8.6435,
                "std": 0.0825,
                "n": 72
            },
            "TLR9": {
                "mean": 3.9546,
                "std": 0.1261,
                "n": 72
            },
            "MYD88": {
                "mean": 9.0896,
                "std": 0.0757,
                "n": 72
            },
            "IRAK4": {
                "mean": 6.1308,
                "std": 0.1777,
                "n": 72
            },
            "TICAM1": {
                "mean": 5.4453,
                "std": 0.0603,
                "n": 72
            },
            "TNFSF13B": {
                "mean": 10.4293,
                "std": 0.1019,
                "n": 72
            },
            "CD40LG": {
                "mean": 6.2525,
                "std": 0.1222,
                "n": 72
            },
            "IRF7": {
                "mean": 7.1697,
                "std": 0.0999,
                "n": 72
            },
            "STAT1": {
                "mean": 10.0714,
                "std": 0.1258,
                "n": 72
            },
            "STAT3": {
                "mean": 8.1525,
                "std": 0.0168,
                "n": 72
            },
            "IFIT1": {
                "mean": 6.7945,
                "std": 0.2825,
                "n": 72
            },
            "IFIT3": {
                "mean": 8.0489,
                "std": 0.2247,
                "n": 72
            },
            "ISG15": {
                "mean": 9.2776,
                "std": 0.2632,
                "n": 72
            },
            "IFI44": {
                "mean": 8.6596,
                "std": 0.2589,
                "n": 72
            },
            "IFI44L": {
                "mean": 7.7985,
                "std": 0.3314,
                "n": 72
            },
            "IFI6": {
                "mean": 9.8696,
                "std": 0.1776,
                "n": 72
            },
            "OAS1": {
                "mean": 7.3119,
                "std": 0.1355,
                "n": 72
            },
            "OAS2": {
                "mean": 5.846,
                "std": 0.2109,
                "n": 72
            },
            "OAS3": {
                "mean": 5.4608,
                "std": 0.1385,
                "n": 72
            },
            "OASL": {
                "mean": 5.9335,
                "std": 0.1138,
                "n": 72
            },
            "RSAD2": {
                "mean": 6.7225,
                "std": 0.3503,
                "n": 72
            },
            "HERC5": {
                "mean": 9.0962,
                "std": 0.2623,
                "n": 72
            },
            "HERC6": {
                "mean": 6.7926,
                "std": 0.1593,
                "n": 72
            },
            "XAF1": {
                "mean": 8.0959,
                "std": 0.2076,
                "n": 72
            },
            "BATF2": {
                "mean": 3.6034,
                "std": 0.1572,
                "n": 72
            },
            "EPSTI1": {
                "mean": 8.9517,
                "std": 0.2287,
                "n": 72
            },
            "DDX58": {
                "mean": 6.2131,
                "std": 0.0477,
                "n": 72
            },
            "DDX60": {
                "mean": 6.8178,
                "std": 0.1448,
                "n": 72
            },
            "IFIH1": {
                "mean": 7.7745,
                "std": 0.105,
                "n": 72
            },
            "CMPK2": {
                "mean": 3.8264,
                "std": 0.1849,
                "n": 72
            },
            "IFITM1": {
                "mean": 11.7854,
                "std": 0.0827,
                "n": 72
            },
            "CXCL10": {
                "mean": 5.4835,
                "std": 0.281,
                "n": 72
            },
            "CCL5": {
                "mean": 12.0922,
                "std": 0.0582,
                "n": 72
            },
            "TNF": {
                "mean": 7.0239,
                "std": 0.0666,
                "n": 72
            },
            "IL1B": {
                "mean": 8.9167,
                "std": 0.1474,
                "n": 72
            },
            "IL6": {
                "mean": 3.3717,
                "std": 0.0175,
                "n": 72
            },
            "ELANE": {
                "mean": 5.1652,
                "std": 0.2137,
                "n": 72
            },
            "ESR1": {
                "mean": 3.3969,
                "std": 0.059,
                "n": 72
            },
            "CYP19A1": {
                "mean": 3.3219,
                "std": 0.0,
                "n": 72
            },
            "IFNA1": {
                "mean": 3.3403,
                "std": 0.0137,
                "n": 72
            }
        }
    },
    "Adult": {
        "n": 25,
        "source": "GSE81622 (adult PBMC, 25 samples)",
        "genes": {
            "TLR7": {
                "mean": 446.0362,
                "std": 112.4957,
                "n": 25
            },
            "TLR8": {
                "mean": 653.0714,
                "std": 501.5561,
                "n": 75
            },
            "TLR9": {
                "mean": 146.0203,
                "std": 11.7445,
                "n": 25
            },
            "MYD88": {
                "mean": 1062.1284,
                "std": 238.9283,
                "n": 25
            },
            "IRAK4": {
                "mean": 326.9169,
                "std": 56.4456,
                "n": 25
            },
            "TICAM1": {
                "mean": 418.6605,
                "std": 272.1604,
                "n": 50
            },
            "TNFSF13B": {
                "mean": 2122.2732,
                "std": 684.8853,
                "n": 50
            },
            "CD40LG": {
                "mean": 279.2183,
                "std": 51.4494,
                "n": 25
            },
            "IRF7": {
                "mean": 738.1696,
                "std": 730.9419,
                "n": 75
            },
            "STAT1": {
                "mean": 2990.7373,
                "std": 1992.7293,
                "n": 75
            },
            "STAT3": {
                "mean": 878.0651,
                "std": 486.0336,
                "n": 75
            },
            "IFIT1": {
                "mean": 525.2465,
                "std": 553.4701,
                "n": 50
            },
            "IFIT3": {
                "mean": 377.6876,
                "std": 253.8374,
                "n": 75
            },
            "ISG15": {
                "mean": 1557.8345,
                "std": 711.5613,
                "n": 25
            },
            "IFI44": {
                "mean": 1782.5948,
                "std": 944.8621,
                "n": 25
            },
            "IFI44L": {
                "mean": 1243.3752,
                "std": 996.7121,
                "n": 25
            },
            "IFI6": {
                "mean": 2190.5921,
                "std": 1545.4701,
                "n": 50
            },
            "MX1": {
                "mean": 4561.5364,
                "std": 1678.8537,
                "n": 25
            },
            "OAS1": {
                "mean": 482.4817,
                "std": 298.5908,
                "n": 100
            },
            "OAS2": {
                "mean": 957.0695,
                "std": 1236.624,
                "n": 100
            },
            "OAS3": {
                "mean": 378.9341,
                "std": 314.3763,
                "n": 50
            },
            "OASL": {
                "mean": 288.0315,
                "std": 131.8865,
                "n": 50
            },
            "RSAD2": {
                "mean": 273.7047,
                "std": 112.4528,
                "n": 25
            },
            "HERC5": {
                "mean": 1606.9617,
                "std": 578.427,
                "n": 25
            },
            "HERC6": {
                "mean": 626.1557,
                "std": 189.5589,
                "n": 25
            },
            "XAF1": {
                "mean": 1246.2323,
                "std": 794.6841,
                "n": 50
            },
            "BATF2": {
                "mean": 138.6241,
                "std": 10.2338,
                "n": 25
            },
            "EPSTI1": {
                "mean": 1668.5847,
                "std": 655.979,
                "n": 25
            },
            "DDX58": {
                "mean": 251.9153,
                "std": 42.2316,
                "n": 25
            },
            "DDX60": {
                "mean": 488.3311,
                "std": 122.3431,
                "n": 25
            },
            "IFIH1": {
                "mean": 563.1244,
                "std": 107.8646,
                "n": 25
            },
            "CMPK2": {
                "mean": 160.8246,
                "std": 29.5355,
                "n": 25
            },
            "IFITM1": {
                "mean": 5441.7478,
                "std": 1179.7176,
                "n": 25
            },
            "CXCL10": {
                "mean": 325.9074,
                "std": 106.8332,
                "n": 25
            },
            "CCL5": {
                "mean": 14076.9991,
                "std": 5977.4652,
                "n": 50
            },
            "TNF": {
                "mean": 1152.7641,
                "std": 1409.2324,
                "n": 25
            },
            "IL1B": {
                "mean": 2402.3821,
                "std": 4890.2757,
                "n": 25
            },
            "IL6": {
                "mean": 144.4391,
                "std": 48.9235,
                "n": 25
            },
            "ELANE": {
                "mean": 569.7907,
                "std": 1093.843,
                "n": 25
            },
            "ESR1": {
                "mean": 135.3705,
                "std": 9.7801,
                "n": 25
            },
            "CYP19A1": {
                "mean": 125.8921,
                "std": 8.846,
                "n": 75
            },
            "IFNA1": {
                "mean": 120.1939,
                "std": 5.4887,
                "n": 25
            }
        }
    }
}

# Alias for backward compatibility
ref = reference


def compare_to_healthy(patient_values, gene, source="All"):
    """Compare patient expression to healthy reference.
    
    Args:
        patient_values: dict of {sample_id: expression_value}
        gene: gene symbol (e.g. "IRF7")
        source: "All", "Pediatric", or "Adult"
    
    Returns:
        dict with fold changes, z-scores, summary
    """
    if source not in reference:
        # Try merging pediatric + adult
        merged = {}
        for s in reference:
            gd = reference[s]["genes"].get(gene)
            if gd:
                merged.setdefault("mean", 0)
                merged["mean"] += gd["mean"] * gd["n"]
                merged.setdefault("n", 0)
                merged["n"] += gd["n"]
                merged["std"] = gd["std"]
        if not merged:
            return None
        merged["mean"] /= merged["n"]
        healthy = merged
    else:
        healthy = reference[source]["genes"].get(gene)
    
    if not healthy:
        return None
    
    vals = np.array(list(patient_values.values()))
    fc = vals / healthy["mean"]
    std = healthy["std"]
    z = (vals - healthy["mean"]) / std if std > 0 else np.zeros_like(vals)
    
    return {
        "gene": gene,
        "source": source,
        "healthy_mean": healthy["mean"],
        "healthy_std": healthy["std"],
        "healthy_n": healthy["n"],
        "patient_mean": float(np.mean(vals)),
        "patient_n": int(len(vals)),
        "fold_change": float(np.mean(fc)),
        "mean_z_score": float(np.mean(z)),
        "samples": {sid: {
            "expression": float(v),
            "fold_change": float(v / healthy["mean"]),
            "z_score": float((v - healthy["mean"]) / std) if std > 0 else 0.0
        } for sid, v in patient_values.items()}
    }


def multi_gene_report(patient_df, source="All", min_n=3):
    """Generate fold-change report comparing patients to healthy for all genes.
    
    Args:
        patient_df: DataFrame with genes as rows, samples as columns
        source: "All", "Pediatric", or "Adult"
        min_n: minimum healthy samples to include
    
    Returns:
        DataFrame sorted by fold change (descending)
    """
    import pandas as pd
    results = []
    
    for gene in patient_df.index:
        # Get healthy reference (merge sources if needed)
        healthy = None
        if source in reference:
            healthy = reference[source]["genes"].get(gene)
        else:
            merged = {}
            for s in reference:
                gd = reference[s]["genes"].get(gene)
                if gd:
                    merged.setdefault("mean", 0)
                    merged["mean"] += gd["mean"] * gd["n"]
                    merged.setdefault("n", 0)
                    merged["n"] += gd["n"]
                    merged["std"] = gd["std"]
            if merged:
                merged["mean"] /= merged["n"]
                healthy = merged
        
        if not healthy or healthy["n"] < min_n:
            continue
        
        vals = patient_df.loc[gene].dropna().values.astype(float)
        if len(vals) < 3:
            continue
        
        fc = float(np.mean(vals / healthy["mean"]))
        std = healthy["std"]
        z = float(np.mean((vals - healthy["mean"]) / std)) if std > 0 else 0.0
        
        results.append({
            "Gene": gene,
            "Healthy_mean": healthy["mean"],
            "Healthy_std": healthy["std"],
            "Patient_mean": float(np.mean(vals)),
            "Fold_Change": round(fc, 4),
            "Z_Score": round(z, 3),
            "Patient_n": len(vals),
            "Healthy_n": healthy["n"]
        })
    
    return pd.DataFrame(results).sort_values("Fold_Change", ascending=False)
