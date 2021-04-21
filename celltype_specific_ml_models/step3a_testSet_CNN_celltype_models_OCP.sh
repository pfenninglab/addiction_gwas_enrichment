#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
##SBATCH --time=0-12
#SBATCH --job-name=cnn_test
#SBATCH --gres=gpu:1
#SBATCH --mem=45G
#SBATCH --array=1-25
#SBATCH --error=logs/cnn_ocp_testSet_%A_%a_output.txt
#SBATCH --output=logs/cnn_ocp_testSet_%A_%a_output.txt


PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
PROJDIR=/home/bnphan/projects/celltype_specific_ml_models
cd $PROJDIR; mkdir -p $PROJDIR/test_performance
source activate tf2

# Get model and affiliated validation seq
LABEL=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $3}' tables/SingleTask_IDR_bestModels_CNN_performance.tsv`
MODEL=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $4}' tables/SingleTask_IDR_bestModels_CNN_performance.tsv`
CELLTYPE=$(echo $LABEL| cut -d '_' -f1-3 )
PREFIX=$(echo $LABEL| cut -d '_' -f1-4 )

# input validation files for 10-fold cross evalutation
TESTPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_testPos.fa
TESTNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_testNeg10x.fa

# score test sequences with default CNN models with default parameters
if [[ ! -f ${PROJDIR}/test_performance/predictions/${PREFIX}/$(basename $MODEL .h5).performance.feather ]]; then
python train_singleTask_CNN_classifier_OCP.py \
--out_dir $PROJDIR/test_performance --model_name $MODEL --mode 'evaluate' \
--valid_fasta_pos $TESTPOSFILE --valid_fasta_neg $TESTNEGFILE
fi

# get scores for false positive rates
if [[ ! -f ${PROJDIR}/predictions/${PREFIX}/$(basename $MODEL .h5).test_negatives.predictions.txt ]]; then
python train_singleTask_CNN_classifier_OCP.py \
--mode 'predict' --model_name $MODEL --prefix $PREFIX --out_dir $PROJDIR \
--predict_fasta $TESTNEGFILE --predict_out "test_negatives"
fi

