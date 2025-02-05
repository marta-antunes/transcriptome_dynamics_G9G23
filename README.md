# Evolutionary dynamics of gene expression during thermal adaptation
Authors: Marta A. Antunes, Marta A. Santos, Mauro Santos, Margarida Matos & Pedro Simões

# - Workflow -
Below we provide a step-by-step workflow of the analyses in the paper. Shaded background represents the code for the analysis or filenames for specific files containing code (provided).

# 1. Overall gene expression analysis
1.1 normalize counts with DESeq2 within galaxy  
input files: count files for all samples (Generation 9 and Genration 23 samples tested in the control environment)
output file: Galaxy340-Normalized_counts.tabular (file with normalized counts)  
all normalizations were done by selecting “output options”>”Output normalised counts” on DESeq2 within Galaxy

1.2 differential expression analysis
input file: Galaxy340-Normalized_counts.tabular
```
run_glmmTMB1binN_on_ExpressionData_RunforGeneration.R
```
output file:DEG_generation_overallAnalysis.csv (table 1 and supplementary table S6)

# 2. Candidate genes at generations 9 and 23 and evolutionary dynamics between generations
2.1. Candidate genes at generation 9
```
run_glmmTMB1binN_on_PTExpressionData_Generation9.R
```

2.2. Candidate genes at generation 23
```
$run_glmmTMB1binN_on_PTExpressionData_Generation23.R
```

2.3. Evolutionary dynamics between generations

|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |count files for PT_G9_C, WPT_G9_C, PT_G23_C and WPT_G23_C samples	|count files for NL_G9_C, WNL_G9_C, NL_G23_C and WNL_G23_C samples|
|Output                          |	Galaxy371-Normalized_counts.tabular	|Galaxy374-Normalized_counts.tabular|

```
run_glmmTMB1binN_on_generation_PT_2models_ratio.R
```

```
run_glmmTMB1binN_on_generation_NL_2models_ratio.R 
```

# 3. Analysis of each generation separately, including data from both historical populations

```
run_glmmTMB1binN_on_SelectionHistory.R
```

|                                | Generation 9       | Generation 23     |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |Galaxy201-Normalized_counts.tabular	|Galaxy213-Normalized_counts.tabular|
|Output                          |	Selection_history_generation9.csv	|Selection_history_generation23.csv|

# 4. Temporal patterns of changes of differential expression of candidate genes

```
scatter_plot_axis_limited_legendInTheBottom.R
```

# 5.	Plateau analysis


# 6. Gene Set Enrichment Analysis (GSEA)	


