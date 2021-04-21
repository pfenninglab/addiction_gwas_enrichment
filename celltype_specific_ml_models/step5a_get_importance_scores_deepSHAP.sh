#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3,gpu
#SBATCH --time=0-12
#SBATCH --job-name=dSHAP
#SBATCH --gres=gpu:1
#SBATCH --mem=120G
#SBATCH --array=1-25
#SBATCH --error=logs/compute_validation_deepSHAP_scores_%A_%a_output.txt
#SBATCH --output=logs/compute_validation_deepSHAP_scores_%A_%a_output.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
PROJDIR=/home/bnphan/projects/celltype_specific_ml_models
cd $PROJDIR; mkdir -p $PROJDIR/shap
source activate tf2

# Get model and affiliated validation seq
LABEL=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $3}' tables/SingleTask_IDR_bestModels_CNN_performance.tsv`
MODEL=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $4}' tables/SingleTask_IDR_bestModels_CNN_performance.tsv`

# validation fasta files 
VALIDPOSFILE=$PROJDIR/FASTA_CV/${LABEL}_validPos.fa
VALIDNEGFILE=$PROJDIR/FASTA_CV/${LABEL}_validNeg10x.fa

#################################################################################
# 1) score validation sequences with default CNN models with default parameters
if [[ ! -f $PROJDIR/shap/${LABEL}.pos2000.imp_SHAP_scores.txt ]]; then
python deepSHAP_importance_scores.py \
--mode evaluate --model_name $MODEL --predict_out ${LABEL} \
--eval_fasta_pos ${VALIDPOSFILE} --eval_fasta_neg ${VALIDNEGFILE}
fi

##########################################################################
# 2) get the SNP predictions sequences with relative to the 10x negatives
EFFECT_SNP='FASTA/top_addiction_snps_allele_effect_501.fa'
LABEL2=top_addiction_snps_effect_allele.${LABEL}
if [[ ! -f $PROJDIR/shap/${LABEL2}.predict.imp_SHAP_scores.txt ]]; then
python deepSHAP_importance_scores.py \
--mode predict --model_name $MODEL --predict_out ${LABEL2} \
--predict_fasta ${EFFECT_SNP} --eval_fasta_neg ${VALIDNEGFILE}
fi

NONEFF_SNP='FASTA/top_addiction_snps_allele_nonEffect_501.fa'
LABEL2=top_addiction_snps_nonEffect_allele.${LABEL}
if [[ ! -f $PROJDIR/shap/${LABEL2}.predict.imp_SHAP_scores.txt ]]; then
python deepSHAP_importance_scores.py \
--mode predict --model_name $MODEL --predict_out ${LABEL2} \
--predict_fasta ${NONEFF_SNP} --eval_fasta_neg ${VALIDNEGFILE}
fi


