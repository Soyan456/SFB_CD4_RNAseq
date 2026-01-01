import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

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
# Filter for our genes of interest
genes_of_interest = ['Notch1', 'Thy1', 'Cd4', 'Cd8', 'Tlr4', 'Nod2']

# Auto-detect gene column (robust to 'Gene Symbol' vs 'GeneSymbol' etc.)
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

df_filtered = df[df[gene_col].astype(str).str.lower().isin([g.lower() for g in genes_of_interest])].copy()

# Define Sample Groups
ctrl_cols = ["B6CTRL1", "B6CTRL2", "B6CTRL3", "B6CTRL4", "B6CTRL5"]
tet_cols = ["TETRAMER_1", "TETRAMER_2", "TETRAMER_3"]

# Verify sample columns exist in the data
ctrl_present = [c for c in ctrl_cols if c in df_filtered.columns]
tet_present = [c for c in tet_cols if c in df_filtered.columns]
if not ctrl_present or not tet_present:
    raise KeyError(f"Sample columns missing. Found ctrl: {ctrl_present}, tet: {tet_present}")

# Mean expression per group (cast to numeric to avoid object dtypes)
df_filtered['mean_CTRL'] = df_filtered[ctrl_present].apply(pd.to_numeric, errors='coerce').mean(axis=1)
df_filtered['mean_TET'] = df_filtered[tet_present].apply(pd.to_numeric, errors='coerce').mean(axis=1)

# Expression detectability check
# A gene is considered detectably expressed if mean expression > 1
# in at least one group
df_filtered['detectable'] = (
    (df_filtered['mean_CTRL'] > 1) | (df_filtered['mean_TET'] > 1)
)
#Apply detectability filter
df_filtered = df_filtered[df_filtered['detectable']].copy()

# Fold Change Calculation
df_filtered['Fold_Change'] = df_filtered['mean_TET'] / df_filtered['mean_CTRL']
df_filtered['log2FC'] = np.log2(df_filtered['Fold_Change'].replace(0, np.nan))

# T-Test Calculation (use .loc to select columns for each row to avoid AttributeError)
p_values = []
for idx in df_filtered.index:
    ctrl_values = df_filtered.loc[idx, ctrl_present].astype(float).values
    tet_values = df_filtered.loc[idx, tet_present].astype(float).values
    t_stat, p_val = ttest_ind(tet_values, ctrl_values, nan_policy='omit', equal_var=False) # Welch's t-test
    p_values.append(p_val)
df_filtered['p_value'] = p_values

# Multiple Testing Correction with Benjamini-Hochberg
if len(df_filtered['p_value']) > 0:
    reject, pvals_corrected, _, _ = multipletests(df_filtered['p_value'], alpha=0.05, method='fdr_bh')
else:
    reject, pvals_corrected = [], []
df_filtered['p_adj'] = pvals_corrected
df_filtered['significant'] = reject

print("Data analysis complete. Here are the results for the genes of interest:")
print(df_filtered[[gene_col, 'mean_CTRL', 'mean_TET', 'Fold_Change', 'log2FC', 'p_value', 'p_adj', 'significant']])

# Save results
df_filtered.to_excel("thymus_gene_analysis_results.xlsx", index=False)

print(f"Loaded: {file_path.name}")

# --- Plotting: grouped bar plot with individual points ---
# Genes and display order
genes_order = ['Cd4', 'Notch1', 'Thy1', 'Tlr4', 'Nod2']

# Friendly group labels
ctrl_label = 'Non-specific CD4 SP (B6CTRL)'
tet_label = 'SFB-specific CD4 SP (TETRAMER)'

# Build long-form DataFrame for plotting
rows = []
for gene in genes_order:
    # find matching rows (case-insensitive)
    match = df_filtered[df_filtered[gene_col].astype(str).str.lower() == gene.lower()]
    if match.empty:
        continue
    # take first match if duplicates
    row = match.iloc[0]
    # collect control sample values
    for samp in ctrl_present:
        try:
            val = float(row[samp])
        except Exception:
            val = np.nan
        rows.append({'gene': gene, 'group': ctrl_label, 'sample': samp, 'value': val})
    # collect tetramer sample values
    for samp in tet_present:
        try:
            val = float(row[samp])
        except Exception:
            val = np.nan
        rows.append({'gene': gene, 'group': tet_label, 'sample': samp, 'value': val})

plot_df = pd.DataFrame(rows)

sns.set(style='whitegrid')
fig, ax = plt.subplots(figsize=(10, 6))

# Palette: control = blue, tetramer = orange
palette = {ctrl_label: '#4c72b0', tet_label: '#dd8452'}

# Bar plot (means)
sns.barplot(data=plot_df, x='gene', y='value', hue='group', order=genes_order,
            estimator=np.nanmean, ci=None, palette=palette, ax=ax)

# Individual data points
sns.stripplot(data=plot_df, x='gene', y='value', hue='group', order=genes_order,
              dodge=True, jitter=True, color='k', size=6, alpha=0.8, ax=ax, linewidth=0.5)

# Remove duplicate legend entries (stripplot adds a second legend)
handles, labels = ax.get_legend_handles_labels()
# keep only first two handles/labels (bars)
if len(handles) >= 2:
    ax.legend(handles[:2], labels[:2], title='Group')

# Annotate significance for CD4 only if adjusted p < 0.05 (and explicitly do NOT mark Thy1)
cd4_row = df_filtered[df_filtered[gene_col].astype(str).str.lower() == 'cd4']
if not cd4_row.empty:
    p_adj_val = cd4_row['p_adj'].values[0] if 'p_adj' in cd4_row.columns else None
    if p_adj_val is not None and not np.isnan(p_adj_val) and p_adj_val < 0.05:
        # place star above the highest datapoint for CD4
        cd4_vals = plot_df[plot_df['gene'] == 'Cd4']['value']
        y = cd4_vals.max() if not cd4_vals.dropna().empty else 0
        y = y + 0.05 * (plot_df['value'].max() - plot_df['value'].min())
        gene_index = genes_order.index('Cd4')
        ax.text(gene_index, y, '*', ha='center', va='bottom', color='red', fontsize=20)

# Titles and labels
ax.set_xlabel('Gene')
ax.set_ylabel('Normalized expression')
ax.set_xticklabels(genes_order)
plt.tight_layout()

# Save figure
fig_path = script_dir / 'thymus_gene_plot.png'
fig.savefig(fig_path, dpi=300)
print(f"Saved plot to: {fig_path}")
print(df_filtered[[gene_col, 'mean_CTRL', 'mean_TET', 'p_adj']])
