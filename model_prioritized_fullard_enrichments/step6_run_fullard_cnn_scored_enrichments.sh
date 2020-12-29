#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --job-name=ldcts
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/run_enrich_%A_%a.err.txt
#SBATCH --output=logs/run_enrich_%A_%a.out.txt
#SBATCH --array=1-18

source activate ldsc

# get the GWAS for this array job
GWAS=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ../mouse_enrichment/gwas_list_sumstats.txt)
GWAS_Label=$(basename $GWAS | sed 's/.sumstats.gz//g')

# run LD score regression on CNN scored peaks
OUTDIR=ldsc_enrichments/Fullard_cnn_scored_enrichments2
mkdir -p $OUTDIR
python /home/bnphan/src/ldsc/ldsc.py \
	--ref-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/baseline_v1.2/baseline. \
	--w-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/weights/weights.hm3_noMHC. \
	--ref-ld-chr-cts Fullard_cnn_scored_peaks.ldcts \
	--h2-cts $GWAS \
	--out $OUTDIR/Fullard_CNNscored_enrichment_${GWAS_Label}

# # run LD score regression on SVM scored peaks
# OUTDIR=ldsc_enrichments/Fullard_svm_scored_enrichments
# python /home/bnphan/src/ldsc/ldsc.py \
# 	--ref-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/baseline_v1.2/baseline. \
# 	--w-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/weights/weights.hm3_noMHC. \
# 	--ref-ld-chr-cts Fullard_svm_scored_peaks.ldcts \
# 	--h2-cts $GWAS \
# 	--out $OUTDIR/Fullard_SVMscored_enrichment_${GWAS_Label}

# # run LD score regression on CNN and SVM combined scored peaks
# OUTDIR=ldsc_enrichments/Fullard_both_scored_enrichments
# python /home/bnphan/src/ldsc/ldsc.py \
# 	--ref-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/baseline_v1.2/baseline. \
# 	--w-ld-chr /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/weights/weights.hm3_noMHC. \
# 	--ref-ld-chr-cts Fullard_both_scored_peaks.ldcts \
# 	--h2-cts $GWAS \
# 	--out $OUTDIR/Fullard_BOTHscored_enrichment_${GWAS_Label}

