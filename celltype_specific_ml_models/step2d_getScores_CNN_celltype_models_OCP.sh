#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
#SBATCH --time=0-12
#SBATCH --job-name=getScores
#SBATCH --gres=gpu:1
#SBATCH --mem=100G
#SBATCH --array=1-55%5
#SBATCH --error=logs/get_scores_cnnv3_%A_%a_output.txt
#SBATCH --output=logs/get_scores_cnnv3_%A_%a_output.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
cd $PROJDIR; source activate tf2

# for SLURM_ARRAY_TASK_ID in {1..55}; do
## cell type and fold index
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $1}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`
FOLD=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $2}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`
echo "Working on ${CELLTYPE}_${FOLD}."

## input training files 
TRAINPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainPos.fa
TRAINNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainNeg10x.fa

## input validation files for 10-fold cross evalutation
VALIDPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validPos.fa
VALIDNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validNeg10x.fa
PREFIX=${CELLTYPE}_${FOLD}

## cyclical learning rate parameters
BATCH_SIZE=1000
EPOCHS=23
NUM_CYCLES=2.35
BASE_LR=1e-2
MAX_LR=1e-1
BASE_M=.85
MAX_M=.99
L2_REG=1e-10
DROPOUT=.25

## score validation sequences with default CNN models with default parameters
if [[ ! -f ${PROJDIR}/predictions/${PREFIX}/${PREFIX}_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.25.validation_positives.predictions.txt ]]; then
python train_singleTask_CNN_classifier_OCP.py --mode 'predict' \
	--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
	--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
	--predict_fasta $VALIDPOSFILE --prefix $PREFIX --out_dir $PROJDIR --predict_out "validation_positives"
fi

if [[ ! -f ${PROJDIR}/predictions/${PREFIX}/${PREFIX}_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.25.validation_negatives.predictions.txt ]]; then
python train_singleTask_CNN_classifier_OCP.py --mode 'predict' \
	--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
	--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
	--prefix $PREFIX --out_dir $PROJDIR \
	--predict_fasta $VALIDNEGFILE --predict_out "validation_negatives"
fi 

