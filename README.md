# Bacterial proliferation in (spinach) leaf material

We estimate the proliferation of bacteria (Pseudomonas) in spinach leaf material by extracting bacteria from leaf discs. 

## Input: TEMPLATE_FILE_CFU_count_calculations.xlsx 
### Explanation: 
- All CFU count scripts use this file as an input
- In this excel sheet we collect all the data of an experiment. E.g. strain, strain code, etc. 
- We do some simple calculations in this sheet to calculat based on the CFU counts. Find further explanations in the "/melanorian/Lab_Protocols_MelanieMendel"

## script CFU_count_calculations_DC3000_library_boxplot.R
### Input:  
TEMPLATE_FILE_CFU_count_calculations.xlsx  type file (CFU_count_calculations_DC3000_library.xlsx)

### Analysis steps:
- log10 transromation of CFU counts
- normalisation of CFU counts based on negative control D36E
- scaling of results based on the positive control DC3000
- statistical analysis after transromation/normalisation/scaling steps using ANOVA + TUKEY
- satistical analysis alternative Dunnet's test per batch using D36E as reference sample
- Outlier removal based on 1.5 x IQR rule 

### Output:
- Plotting of results at all stages in a boxplot + scatter 
- Plotting possible with/without statistics
- Summary .txt file 
- possibility to export data after each step
- possiblity to export results as input for plotting in a summary heat map

## script CFU_count_calculations_grouped_boxplots_tt.R
### Input:  
TEMPLATE_FILE_CFU_count_calculations.xlsx  type file (CFU_count_calculations_Pseudomonas_strains.xlsx)

### Analysis steps:
- log10 transromation of CFU counts
- no further normalisation
- Statistical analysis using a two sided t-test comparing x dpi and y dpi 

### Output:
- Plotting of results as a grouped boxplot + scatter
- each boxplot showes the results for one bacterial strain
- each boxplot showes the results for 3 different ODs (0.2, 0.4, 0.6)
- each boxplot showes the results for x dpi and y dpi (can be defined, here 0 dpi vs 2 dpi
- boxplots are grouped by dpi
- Statistical results are safed in PreDate_tt_summary_CUF_Pst_strains.txt
