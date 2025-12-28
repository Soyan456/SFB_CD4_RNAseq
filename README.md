# SFB_CD4_RNAseq=
RNA-seq Analysis of SFB-Specific vs Non-Specific Thymic CD4 Single-Positive T Cells

**This repository contains code for analyzing publicly available RNA-seq data to compare gene expression in SFB-specific versus non-specific thymic CD4 single-positive (SP) T cells in Mus musculus. The analysis focuses on key genes involved in thymic T cell development (Cd4, Cd8, Notch1, Thy1) and microbiota recognition pathways (Tlr4, Nod2). Differential expression between groups was assessed using Welch’s t-tests with Benjamini-Hochberg correction for multiple testing.**

**Data Source**: The RNA-seq dataset was obtained from the Gene Expression Omnibus (GEO) database (accession number: GSE171279). The dataset includes normalized counts for:
1. SFB-tetramer–positive thymic CD4 T cells
2. Non-specific (control) thymic CD4 T cells
(Each sample represents pooled thymic tissue from ten mice.)

**Requirements / Dependencies**: The code was written in Python 3.12.10 and uses the following libraries: pandas, numpy, scipy, statsmodels, matplotlib, seaborn, pathlib
Install dependencies using pip if needed: pip3 install pandas numpy scipy statsmodels matplotlib seaborn

**Files**

"thymus_cd4_rnaseq_analysis.py" — Main script performing the differential expression analysis and generating plots.
"thymus_gene_analysis_results.xlsx" — Output Excel file containing ID, gene symbol, B6CTRL#, TETRAMER#, mean expression, detectability, fold change, log2 fold change, raw p-values, adjusted p-values, and significance.

"thymus_gene_plot.png" — Generated figure showing normalized expression values with individual biological replicates.

**Usage**

Place the RNA-seq data file (CSV or Excel) in the same folder as the script. The script automatically detects the first CSV/XLS/XLSX file.

Run the script: python thymus_cd4_rnaseq_analysis.py

The results table and bar plot will be saved in the same folder as the script.

**Analysis Workflow**
1. Filter genes of interest: Cd4, Cd8, Notch1, Thy1, Tlr4, Nod2
2. Calculate mean expression for control and SFB-specific groups
3. Determine detectability (mean expression > 1)
4. Calculate fold change and log2 fold change
5. Perform Welch’s two-sample t-test
6. Adjust p-values using Benjamini-Hochberg correction
7. Generate a grouped bar plot with overlaid biological replicates
8. Annotate statistically significant differences for Cd4

**Results**
Cd4 expression is significantly lower in SFB-specific thymic CD4 SP T cells compared to controls
Notch1, Thy1, Tlr4, and Nod2 show no statistically significant differences
Fold change, log2 fold change, and adjusted p-values are provided in the results Excel file
The bar plot visualizes mean expression with individual data points

**References**
Gene roles in thymic T cell development: Cd4, Cd8, Notch1, Thy1
Microbiota recognition pathways: Tlr4, Nod2
GEO dataset: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171279]
Python libraries: pandas, numpy, scipy, statsmodels, matplotlib, seaborn

**Notes**
Ensure the sample columns in your data match the names defined in the script (B6CTRL1-5, TETRAMER_1-3). Adjust if necessary.
Only the first Excel/CSV file in the folder will be processed.
CD8 gene is not included in the dataset, so it is not analyzed in this script.
