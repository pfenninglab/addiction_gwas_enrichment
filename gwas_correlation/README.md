# Comparison of genetic architecture of traits defined by GWAS 

## Dependencies

The following tools are necessary in order to reproduce these analyses. 
- Python 3
  - [ldsc](https://github.com/bulik/ldsc)
- R version 3.6.3
  - ggplot2
  - ggcorrplot
  - tidyverse
  - data.table
  - reshape2
- Slurm Workload Manager

We note that while ```step1_run_gwas_correlation.sh``` is configured to submit a batch job using Slurm, the script can be reconfigured to run each analysis sequentially.

