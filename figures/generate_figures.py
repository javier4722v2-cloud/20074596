#!/usr/bin/env python3
"""
Generate figures for the lupus paper.
Output: lumus/figures/fig1_pathway.png, fig2_heatmap.png, fig3_blood_vs_skin.png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os

OUT = r'C:\Users\javie\.openclaw\workspace\lumus\figures'
os.makedirs(OUT, exist_ok=True)

# ============================================================
# COLOUR PALETTE (colourblind-friendly)
# ============================================================
BLOOD_COL = '#2166AC'   # blue - blood
SKIN_COL  = '#D6604D'    # red - skin
CTRL_COL  = '#878787'    # grey - control
APS_COL   = '#4DAF4A'    # green - APS (via WGCNA)
GOLD_COL  = '#FFD700'    # gold - healthy

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ============================================================
# FIGURE 1: Conceptual Pathway Diagram
# ============================================================
def fig1_pathway():
    fig, ax = plt.subplots(1, 1, figsize=(12, 7))
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1, 8)
    ax.axis('off')
    ax.set_title('Figure 1. The TLR7\u2192BAFF Self-Sustaining Loop:\nBlood (Systemic) vs Skin (Local) IFN Activation',
                 fontsize=14, fontweight='bold', pad=10)

    # Helper: draw a box with text
    def box(x, y, text, w=2.0, h=0.8, face='#E6F0FA', edge='#2166AC', fontsize=10, bold=False):
        r = mpatches.FancyBboxPatch((x-w/2, y-h/2), w, h,
                                     boxstyle="round,pad=0.1",
                                     facecolor=face, edgecolor=edge, linewidth=2)
        ax.add_patch(r)
        fw = 'bold' if bold else 'normal'
        ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
                fontweight=fw, color='#333333')

    def arrow(x1, y1, x2, y2, label='', col='#555555', ls='-', lw=1.8):
        dx, dy = x2-x1, y2-y1
        ax.annotate('', xy=(x2-0.08*dx, y2-0.08*dy), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color=col, lw=lw,
                                    linestyle=ls, connectionstyle='arc3,rad=0'))
        if label:
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax.text(mx, my+0.05, label, ha='center', va='bottom',
                    fontsize=8, style='italic', color=col)

    # === LEFT: BLOOD LOOP ===
    ax.text(1.5, 7.5, 'BLOOD (all lupus variants)', fontsize=12, fontweight='bold',
            color=BLOOD_COL, ha='center')

    # TLR7 (top left)
    box(1, 6.2, 'TLR7', w=1.3, h=0.6, face='#D4E6F1', edge=BLOOD_COL, bold=True)
    box(1, 5.0, 'MYD88', w=1.3, h=0.6, face='#D4E6F1', edge=BLOOD_COL)
    box(1, 3.8, 'IRF7', w=1.3, h=0.6, face='#A9CCE3', edge=BLOOD_COL, bold=True)

    arrow(1, 5.9, 1, 5.3)

    # CXCL10 + IFN sig downstream
    box(1, 2.2, 'IFN signature\n(IFIT1, ISG15,\nIFI44, CXCL10)', w=2.2, h=1.2,
        face='#E8F8F5', edge=BLOOD_COL)
    arrow(1, 3.5, 1, 2.8)

    # BAFF branch
    box(4, 5.5, 'TNFSF13B\n(BAFF)', w=1.5, h=0.8, face='#D4E6F1', edge='#1B4F72')
    box(4, 4.0, 'B-cell\nsurvival', w=1.3, h=0.7, face='#EBF5FB', edge='#1B4F72')
    arrow(2.3, 5.5, 3.3, 5.5, 'NF-kB')
    arrow(4, 4.75, 4, 5.1)

    # BAFF -> IRF7 feedback
    arrow(4, 3.1, 1.65, 3.5, 'BAFF \u2192 IRF7', col=BLOOD_COL, ls='--')

    # IFN -> TLR7 feed-forward
    arrow(1, 1.6, 1, 0.5, '', col=BLOOD_COL)
    box(1, -0.1, 'Autoantibodies\n\u2191 RNA immune complexes', w=2.2, h=0.8,
        face='#F9E79F', edge='#B7950B')
    arrow(1, 0.3, 1, 0.7, 'loop closes', col='#B7950B', ls='--')

    # Autoantibodies back up to TLR7
    ax.annotate('', xy=(1.65, 5.9), xytext=(0.35, 0.4),
                arrowprops=dict(arrowstyle='->', color='#B7950B', lw=1.2,
                                linestyle='dotted', connectionstyle='arc3,rad=0.3'))
    ax.text(1.2, 3.0, 'Self-sustaining', fontsize=9, fontweight='bold',
            color=BLOOD_COL, rotation=90, alpha=0.5)

    # Correlation numbers
    ax.text(6, 4.5, 'Correlations (n=924 SLE):', fontsize=9, fontweight='bold', color='#444')
    ax.text(6, 4.1, 'IRF7\u2192IFIT1:  r = +0.86', fontsize=9, color=BLOOD_COL)
    ax.text(6, 3.8, 'IRF7\u2192ISG15: r = +0.91', fontsize=9, color=BLOOD_COL)
    ax.text(6, 3.5, 'BAFF\u2192IRF7:  r = +0.74', fontsize=9, color=BLOOD_COL)
    ax.text(6, 3.2, 'TLR7\u2192IRF7:  r = +0.25', fontsize=9, color='#888')
    ax.text(6, 2.9, '(IRF7 receives multiple inputs)', fontsize=8, style='italic', color='#888')

    # === RIGHT: SKIN MECHANISM ===
    ax.text(9, 7.5, 'SKIN (CLE lesions)', fontsize=12, fontweight='bold',
            color=SKIN_COL, ha='center')

    box(9, 6.0, 'UV / Damage', w=1.3, h=0.6, face='#FADBD8', edge=SKIN_COL, bold=True)
    box(9, 4.6, 'cGAS-STING\n(alternative IFN)', w=1.8, h=0.9, face='#FADBD8', edge=SKIN_COL)
    box(9, 3.2, 'IFN-\u03ba\n(keratinocyte-specific)', w=1.8, h=0.8,
        face='#F5B7B1', edge=SKIN_COL, bold=True)
    box(9, 1.8, 'IFN signature\n(CXCL10, ISG15,\nIFI44, etc.)', w=2.0, h=1.0,
        face='#FDEDEC', edge=SKIN_COL)

    arrow(9, 5.7, 9, 5.1)
    arrow(9, 4.15, 9, 3.6)
    arrow(9, 2.8, 9, 2.3)

    # X mark over TLR7 in skin
    ax.text(7.5, 6.2, '\u2718 TLR7/BAFF disconnected',
            fontsize=10, color=SKIN_COL, fontweight='bold')
    ax.text(7.5, 5.85, 'TLR7\u2192BAFF: r = -0.18',
            fontsize=9, color='#888')
    ax.text(7.5, 5.6, 'IRF7\u2192CXCL10: r = +0.12',
            fontsize=9, color='#888')

    fig.tight_layout()
    path = os.path.join(OUT, 'fig1_pathway.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f'Saved: {path}')
    return path


# ============================================================
# FIGURE 2: Multi-Type Correlation Heatmap
# ============================================================
def fig2_heatmap():
    # Load mega correlation data
    path = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data\results\mega_correlation_final.csv'
    df = pd.read_csv(path)

    # Select groups to display (blood lupus types + CLE skin + controls)
    groups_of_interest = [
        'SLE blood (GSE65391)',
        'SLE blood (GSE121239)',
        'SLE blood (GSE81622)',
        'CLE blood',
        'CLE skin (all)',
        'Normal skin',
        'Psoriasis skin',
        'Control blood (GSE65391)',
        'Control blood (GSE81622)',
    ]

    # Define display labels
    display_labels = [
        'SLE paediatric\nblood (n=924)',
        'SLE adult\nblood (n=312)',
        'SLE adult\nblood (n=28)',
        'CLE\nblood (n=62)',
        'CLE\nskin (n=25)',
        'Normal\nskin (n=14)',
        'Psoriasis\nskin (n=17)',
        'Healthy\nblood (n=72)',
        'Healthy\nblood (n=27)',
    ]

    gene_pairs = ['IRF7->IFIT1', 'IRF7->ISG15', 'TNFSF13B->IRF7',
                  'IRF7->CXCL10', 'TLR7->CXCL10', 'MYD88->TNFSF13B',
                  'TLR7->IRF7', 'TLR7->TNFSF13B', 'ELANE->CXCL10']

    pair_labels = ['IRF7\u2192IFIT1', 'IRF7\u2192ISG15', 'BAFF\u2192IRF7',
                   'IRF7\u2192CXCL10', 'TLR7\u2192CXCL10', 'MYD88\u2192BAFF',
                   'TLR7\u2192IRF7', 'TLR7\u2192BAFF', 'ELANE\u2192CXCL10']

    # Build matrix
    n_groups = len(groups_of_interest)
    n_pairs = len(gene_pairs)
    r_mat = np.full((n_groups, n_pairs), np.nan)

    for i, g in enumerate(groups_of_interest):
        row = df[df['group'] == g]
        if len(row) == 0:
            continue
        for j, pair in enumerate(gene_pairs):
            r_col = f'{pair}_r'
            if r_col in row.columns:
                r_mat[i, j] = row.iloc[0][r_col]

    # Plot
    fig, ax = plt.subplots(figsize=(11, 6.5))
    cmap = plt.cm.RdBu_r

    # Create masked array for display
    masked = np.ma.masked_invalid(r_mat)
    im = ax.imshow(masked, cmap=cmap, vmin=-0.5, vmax=1.0, aspect='auto')

    # Add text annotations
    for i in range(n_groups):
        for j in range(n_pairs):
            val = r_mat[i, j]
            if not np.isnan(val):
                txt = f'{val:.2f}'
                col = 'white' if abs(val) > 0.55 else '#333'
                ax.text(j, i, txt, ha='center', va='center', fontsize=7, color=col)

    # Labels
    ax.set_xticks(range(n_pairs))
    ax.set_xticklabels(pair_labels, rotation=30, ha='right', fontsize=9)
    ax.set_yticks(range(n_groups))
    ax.set_yticklabels(display_labels, fontsize=9)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('Pearson r', fontsize=10)

    # Separate blood from skin/controls with lines
    # (groups 0-3 are blood lupus, 4-6 are skin, 7-8 are healthy blood)
    ax.axhline(y=3.5, color='#444', linewidth=2, linestyle='-')
    ax.axhline(y=6.5, color='#444', linewidth=2, linestyle='-')

    # Section labels
    ax.text(-0.6, 1.7, 'LUPUS\nBLOOD', fontsize=10, fontweight='bold', color=BLOOD_COL,
            ha='center', va='center', rotation=90)
    ax.text(-0.6, 5.2, 'SKIN\n(control)', fontsize=10, fontweight='bold', color=SKIN_COL,
            ha='center', va='center', rotation=90)
    # Also add psoriasis in skin
    ax.text(-0.6, 7.2, 'HEALTHY\nBLOOD', fontsize=10, fontweight='bold', color=CTRL_COL,
            ha='center', va='center', rotation=90)

    ax.set_title('Figure 2. Loop Correlation Strength Across Tissue Types and Disease States',
                 fontsize=13, fontweight='bold', pad=8)
    fig.tight_layout()
    path = os.path.join(OUT, 'fig2_heatmap.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f'Saved: {path}')
    return path


# ============================================================
# FIGURE 3: Blood vs Skin Correlation Comparison
# ============================================================
def fig3_blood_vs_skin():
    path = r'C:\Users\javie\.openclaw\workspace\lumus\cle_data\results\mega_correlation_final.csv'
    df = pd.read_csv(path)

    def get_r(group_name, pair):
        row = df[df['group'] == group_name]
        if len(row) == 0:
            return np.nan
        col = f'{pair}_r'
        return row.iloc[0][col] if col in row.columns else np.nan

    # Groups for blood lupus
    blood_lupus = ['SLE blood (GSE65391)', 'SLE blood (GSE121239)',
                   'SLE blood (GSE81622)', 'CLE blood']
    blood_labels = ['SLE ped\n(GSE65391)', 'SLE adult\n(GSE121239)',
                    'SLE adult\n(GSE81622)', 'CLE\nblood']

    # Skin groups
    skin_groups = ['CLE skin (all)', 'Normal skin', 'Psoriasis skin']
    skin_labels = ['CLE\nskin', 'Normal\nskin', 'Psoriasis\nskin']

    pairs = ['IRF7->IFIT1', 'TNFSF13B->IRF7', 'IRF7->CXCL10', 'TLR7->CXCL10']
    pair_labels = ['IRF7\u2192IFIT1', 'BAFF\u2192IRF7', 'IRF7\u2192CXCL10', 'TLR7\u2192CXCL10']

    n_groups = len(blood_lupus) + len(skin_groups)
    n_pairs = len(pairs)

    x = np.arange(n_groups)
    width = 0.18
    multiplier = 0

    fig, ax = plt.subplots(figsize=(11, 5.5))

    for i, (pair, plabel) in enumerate(zip(pairs, pair_labels)):
        offset = width * multiplier
        vals = [get_r(g, pair) for g in blood_lupus + skin_groups]
        color = plt.cm.viridis(i / n_pairs)
        bars = ax.bar(x + offset, vals, width, label=plabel, color=color, edgecolor='#333',
                      linewidth=0.5)
        multiplier += 1

    # Labels
    all_labels = blood_labels + skin_labels
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(all_labels, fontsize=9)
    ax.set_ylabel('Pearson correlation (r)', fontsize=11)
    ax.set_ylim(-0.4, 1.05)
    ax.axhline(y=0, color='#888', linewidth=0.8)

    # Separator between blood and skin
    sep_x = len(blood_lupus) - 0.5 + width * 2
    ax.axvline(x=sep_x, color='#444', linewidth=2, linestyle='--')
    ax.text(-0.3, 1.02, 'LUPUS BLOOD', fontsize=11, fontweight='bold', color=BLOOD_COL)
    ax.text(sep_x + 0.4, 1.02, 'SKIN', fontsize=11, fontweight='bold', color=SKIN_COL)

    # Add significance notes for skin
    ax.annotate('TLR7/BAFF\ndisconnected\nin skin',
                xy=(6.8, 0.1), fontsize=8, color=SKIN_COL, fontweight='bold',
                ha='center')

    ax.legend(loc='lower left', fontsize=9, ncol=2, framealpha=0.9)
    ax.set_title('Figure 3. Dichotomy of Loop Connectivity: Blood vs Skin',
                 fontsize=13, fontweight='bold', pad=8)
    fig.tight_layout()
    path = os.path.join(OUT, 'fig3_blood_vs_skin.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f'Saved: {path}')
    return path


# ============================================================
# FIGURE 4: Organ Invariance Bar Chart (from paper text)
# ============================================================
def fig4_organ_invariance():
    """Show that loop correlations are the same regardless of organ affected."""
    # We'll plot a conceptual bar chart with the data from the paper
    organs = ['All SLE\n(n=924)', 'Renal\n(n=408)', 'Musc\n(n=72)',
              'Mucocut\n(n=67)', 'Haematol\n(n=45)', 'CNS\n(n=7)']
    pairs = {'IRF7\u2192IFIT1': [0.86, 0.86, 0.87, 0.90, 0.88, 0.85],
             'IRF7\u2192ISG15': [0.91, 0.91, 0.89, 0.93, 0.87, 0.98],
             'BAFF\u2192IRF7': [0.77, 0.79, 0.75, 0.68, 0.71, 0.93]}

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(organs))
    width = 0.22
    colors = ['#2166AC', '#4393C3', '#92C5DE']

    for i, (label, vals) in enumerate(pairs.items()):
        offset = (i - 1) * width
        bars = ax.bar(x + offset, vals, width, label=label,
                      color=colors[i], edgecolor='#333', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(organs, fontsize=9)
    ax.set_ylabel('Pearson r', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc='lower right')
    ax.set_title('Figure 4. Loop Correlations Are Invariant Across Organ Involvement (SLEDAI)',
                 fontsize=12, fontweight='bold', pad=8)

    # Add annotation about invariance
    ax.text(0.5, -0.15,
            'Fold-change difference between "affected" and "unaffected" organs < 15% for all genes',
            transform=ax.transAxes, fontsize=9, style='italic', ha='center', color='#555')

    fig.tight_layout()
    path = os.path.join(OUT, 'fig4_organ_invariance.png')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f'Saved: {path}')
    return path


# ============================================================
# Execute all
# ============================================================
if __name__ == '__main__':
    fig1_pathway()
    fig2_heatmap()
    fig3_blood_vs_skin()
    fig4_organ_invariance()
    print('\nAll figures generated successfully.')
