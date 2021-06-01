## Comparison of genetic architecture of traits defined by GWAS 

### Dependencies

The following packages are necessary in order to reproduce these analyses. 
- (ldsc)[https://github.com/bulik/ldsc] 
- ggplot2
- ggcorrplot
- tidyverse
- data.table
- reshape2

```step1_run_gwas_correlation.sh``` is configured to submit a batch job using Slurm Workload Manager, but can be reconfigured to run each analysis sequentially.

```step2_plot_gwas_correlation.R``` is configured to run using R version 3.6.3.