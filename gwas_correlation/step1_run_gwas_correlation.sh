#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --job-name=gwasCorr
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/run_correlate_GWAS_%A_%a.txt
#SBATCH --output=logs/run_correlate_GWAS_%A_%a.txt
#SBATCH --array=1-18

source activate ldsc

# get the GWAS for this array job
GWAS=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/gwas_list_sumstats.txt)
GWAS_Label=$(basename $GWAS | sed 's/.sumstats.gz//g')

# get the other 17 gwas to compute gwas correlation against
OTHER_GWAS=$(awk -v var=${SLURM_ARRAY_TASK_ID} 'BEGIN { ORS = ","; NR!=var}1' /projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/gwas_list_sumstats.txt)
OTHER_GWAS=$(echo ${OTHER_GWAS} | sed 's/\(.*\),/\1 /')

# GWAS correlation w/ all of ith GWAS with all other 17 GWASs
ldsc.py \
--rg ${GWAS},${OTHER_GWAS} \
--ref-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/eur_w_ld_chr/ \
--w-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/eur_w_ld_chr/ \
--out tables/p1_${GWAS_Label}_correlation

