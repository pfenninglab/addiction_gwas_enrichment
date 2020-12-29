#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu
#SBATCH --time=0-12
#SBATCH --job-name=cnn
#SBATCH --gres=gpu:1
#SBATCH --mem=90G
#SBATCH --array=1-55%5
#SBATCH --error=logs/cnn_clr_%A_%a_output.txt
#SBATCH --output=logs/cnn_clr_%A_%a_output.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
cd $PROJDIR

# cell type and fold index
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $1}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`
FOLD=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $2}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`

# input training files 
TRAINPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainPos.fa
TRAINNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainNeg10x.fa

# input validation files for 10-fold cross evalutation
VALIDPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validPos.fa
VALIDNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validNeg10x.fa

LABEL=${CELLTYPE}_${FOLD}
mkdir -p $PROJDIR/ml_cnn/$LABEL

# cyclical learning rate parameters
CLR_MODE='triangular2'
STEP_SIZE_SCALE=4
NUM_CYCLES=8
BASE_LR=2e-05
MAX_LR=.2
L2_REG=1e-10

source activate tf2
# conduct range tests for cyclical learning rate
python train_singleTask_CNN_classifier.py \
	--conv_width 11 --l2_reg $L2_REG \
	--base_lr 1e-6 \
	--max_lr 10 \
	--range_test \
	$LABEL $TRAINPOSFILE $TRAINNEGFILE 

# train CNN models with default parameters
python train_singleTask_CNN_classifier.py \
	--conv_width 11 --l2_reg $L2_REG \
	--clr_mode $CLR_MODE \
	--numCycles $NUM_CYCLES \
	--step_size_scale $STEP_SIZE_SCALE \
	--base_lr $BASE_LR \
	--max_lr $MAX_LR \
	$LABEL $TRAINPOSFILE $TRAINNEGFILE

# score test sequences with default CNN models with default parameters
python train_singleTask_CNN_classifier.py \
	--conv_width 11 --l2_reg $L2_REG \
	--clr_mode $CLR_MODE \
	--numCycles $NUM_CYCLES \
	--step_size_scale $STEP_SIZE_SCALE \
	--base_lr $BASE_LR \
	--max_lr $MAX_LR \
	--predict_mode \
	$LABEL $VALIDPOSFILE $VALIDNEGFILE 






