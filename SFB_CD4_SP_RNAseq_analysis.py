import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# Load the first Excel/CSV file found in the same folder as this script
script_dir = Path(__file__).resolve().parent
candidates = list(script_dir.glob('*.xlsx')) + list(script_dir.glob('*.xls')) + list(script_dir.glob('*.csv'))
if not candidates:
    raise FileNotFoundError(f"No Excel/CSV file found in {script_dir}")

file_path = candidates[0]
if file_path.suffix.lower() == '.csv':
    df = pd.read_csv(file_path, header=1)
else:
    df = pd.read_excel(file_path, header=1)

# Auto-detect gene column
possible_names = {'gene symbol', 'gene_symbol', 'genesymbol', 'gene', 'gene_name', 'symbol', 'gene name'}
gene_col = None
for col in df.columns:
    if col.lower().strip() in possible_names:
        gene_col = col
        break
if gene_col is None:
    for col in df.columns:
        if 'gene' in col.lower():
            gene_col = col
            break
if gene_col is None:
    raise KeyError(f"Could not find a gene column in {file_path.name}. Columns: {list(df.columns)}")

df_filtered = df.copy()

# Normalize gene column immediately to avoid whitespace issues
df_filtered[gene_col] = df_filtered[gene_col].astype(str).str.strip()

# Define Sample Groups
ctrl_cols = ["B6CTRL1", "B6CTRL2", "B6CTRL3", "B6CTRL4", "B6CTRL5"]
tet_cols  = ["TETRAMER_1", "TETRAMER_2", "TETRAMER_3"]

ctrl_present = [c for c in ctrl_cols if c in df_filtered.columns]
tet_present  = [c for c in tet_cols  if c in df_filtered.columns]
if not ctrl_present or not tet_present:
    raise KeyError(f"Sample columns missing. Found ctrl: {ctrl_present}, tet: {tet_present}")

# Mean expression per group
df_filtered['mean_CTRL'] = df_filtered[ctrl_present].apply(pd.to_numeric, errors='coerce').mean(axis=1)
df_filtered['mean_TET']  = df_filtered[tet_present].apply(pd.to_numeric, errors='coerce').mean(axis=1)

# Detectability filter
df_filtered['detectable'] = (df_filtered['mean_CTRL'] > 1) | (df_filtered['mean_TET'] > 1)
df_filtered = df_filtered[df_filtered['detectable']].copy()

# Fold Change
df_filtered['Fold_Change'] = df_filtered['mean_TET'] / df_filtered['mean_CTRL']
df_filtered['log2FC'] = np.log2(df_filtered['Fold_Change'].replace(0, np.nan))

# T-test and BH correction on genes of interest only
genes_of_interest = ['Cd4', 'Notch1', 'Thy1', 'Tlr4', 'Nod2', 'Gapdh']

mask_goi = df_filtered[gene_col].str.lower().isin([g.lower() for g in genes_of_interest])
df_filtered = df_filtered[mask_goi].copy()

p_values = []
for idx in df_filtered.index:
    ctrl_values = df_filtered.loc[idx, ctrl_present].astype(float).values
    tet_values  = df_filtered.loc[idx, tet_present].astype(float).values
    _, p_val = ttest_ind(tet_values, ctrl_values, nan_policy='omit', equal_var=False)
    p_values.append(p_val)
df_filtered['p_value'] = p_values

reject, pvals_corrected, _, _ = multipletests(df_filtered['p_value'], alpha=0.05, method='fdr_bh')
df_filtered['p_adj']       = pvals_corrected
df_filtered['significant'] = reject

# Print results for genes of interest only
mask = df_filtered[gene_col].str.lower().isin([g.lower() for g in genes_of_interest])
print("Data analysis complete. Results for genes of interest:")
print(df_filtered.loc[mask, [gene_col, 'mean_CTRL', 'mean_TET', 'Fold_Change', 'log2FC', 'p_value', 'p_adj', 'significant']])

# Save full results
df_filtered.to_excel("thymus_gene_analysis_results.xlsx", index=False)
print(f"Loaded: {file_path.name}")

# ── PLOTTING ────────────────────────────────────────────────────────────────

genes_order = ['Cd4', 'Notch1', 'Thy1', 'Tlr4', 'Nod2', 'Gapdh']
ctrl_label  = 'Non-specific CD4 SP (B6CTRL)'
tet_label   = 'SFB-specific CD4 SP (TETRAMER)'
palette = {ctrl_label: '#89ABD9', tet_label: '#DD8452'}

# Build plot dataframe (ctrl=1 baseline, tet=fold change)
rows = []
for gene in genes_order:
    match = df_filtered[df_filtered[gene_col].str.lower() == gene.lower()]
    if match.empty:
        continue
    row = match.iloc[0]
    fc  = row['Fold_Change']
    if pd.isna(fc) or np.isinf(fc):
        continue
    rows.append({'gene': gene, 'group': ctrl_label, 'value': 1})
    rows.append({'gene': gene, 'group': tet_label,  'value': fc})

plot_df = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(10, 6))
x     = np.arange(len(genes_order))
width = 0.35

# Draw bars
for group, offset in {ctrl_label: -width/2, tet_label: width/2}.items():
    group_df = plot_df[plot_df['group'] == group]
    values = [
        group_df.loc[group_df['gene'] == g, 'value'].values[0]
        if not group_df[group_df['gene'] == g].empty else np.nan
        for g in genes_order
    ]
    ax.bar(x + offset, values, width, label=group, color=palette[group], alpha=0.85, zorder=2)

# Draw individual sample points
for i, gene in enumerate(genes_order):
    match = df_filtered[df_filtered[gene_col].str.lower() == gene.lower()]
    if match.empty:
        continue
    row = match.iloc[0]
    fc  = row['Fold_Change']
    if pd.isna(fc) or np.isinf(fc):
        continue

    ctrl_vals      = row[ctrl_present].astype(float).values
    ctrl_normed    = np.clip(ctrl_vals / row['mean_CTRL'], 0, 4)
    ax.scatter(
        np.full(len(ctrl_normed), i - width/2) + np.random.uniform(-0.05, 0.05, len(ctrl_normed)),
        ctrl_normed, color='black', s=20, zorder=3, alpha=0.7
    )

    tet_vals    = row[tet_present].astype(float).values
    tet_normed  = np.clip(tet_vals / row['mean_CTRL'], 0, 4)
    ax.scatter(
        np.full(len(tet_normed), i + width/2) + np.random.uniform(-0.05, 0.05, len(tet_normed)),
        tet_normed, color='black', s=20, zorder=3, alpha=0.7
    )

# Significance annotations
sig_genes = {}
for i, gene in enumerate(genes_order):
    match = df_filtered[df_filtered[gene_col].str.lower() == gene.lower()]
    if not match.empty and match.iloc[0]['significant']:
        sig_genes[gene] = (i, match.iloc[0])

plt.draw()
y_max = ax.get_ylim()[1]

for gene, (i, row) in sig_genes.items():
    p_adj = row['p_adj']
    stars = '***' if p_adj < 0.001 else '**' if p_adj < 0.01 else '*'
    
    ctrl_x   = x[i] - width/2   # center of blue bar
    tet_x    = x[i] + width/2   # center of orange bar
    bar_top  = max(1, row['Fold_Change'])
    bracket_y = 1.22
    
    # Draw horizontal bracket line spanning both bars
    ax.plot([ctrl_x, tet_x], [bracket_y, bracket_y], color='black', linewidth=1.2)
    # Short vertical ticks at each end
    ax.plot([ctrl_x, ctrl_x], [bracket_y - 0.03, bracket_y], color='black', linewidth=1.2)
    ax.plot([tet_x,  tet_x],  [bracket_y - 0.03, bracket_y], color='black', linewidth=1.2)
    # Stars centered above the bracket midpoint
    ax.text(
        (ctrl_x + tet_x) / 2, bracket_y + 0.02, stars,
        ha='center', va='bottom',
        fontsize=20, fontweight='bold', color='black'
    )
    
# Axis labels — italicised, first letter uppercase
italic_labels = [rf'$\it{{{g[0].upper() + g[1:]}}}$' for g in genes_order]
ax.set_xticks(x)
ax.set_xticklabels(italic_labels, fontsize=14)
ax.set_xlabel('Gene', fontsize=16)
ax.set_ylabel('Fold Change vs Control', fontsize=16)
ax.set_title('Thymic Gene Expression: SFB-Specific vs Non-Specific CD4 SP T Cells', fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.grid(False)
ax.legend(fontsize=11, loc='upper left', frameon=False)

plt.tight_layout()
fig_path = script_dir / 'thymus_gene_plot.png'
fig.savefig(fig_path, dpi=300)
print(f"Saved plot to: {fig_path}")
