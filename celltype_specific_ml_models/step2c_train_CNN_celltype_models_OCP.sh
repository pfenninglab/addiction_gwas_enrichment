#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu,pfen3
#SBATCH --time=0-12
#SBATCH --job-name=cnn_train
#SBATCH --gres=gpu:1
#SBATCH --mem=80G
#SBATCH --array=1-55%5
#SBATCH --error=logs/train_cnnv3_%A_%a_output.txt
#SBATCH --output=logs/train_cnnv3_%A_%a_output.txt

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
DROPOUT=.2

## train CNN models with default parameters
for DROPOUT in .2 .25; do
python train_singleTask_CNN_classifier_OCP.py --mode 'train' --out_dir $PROJDIR \
	--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
	--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
	--train_fasta_pos $TRAINPOSFILE --train_fasta_neg $TRAINNEGFILE \
	--valid_fasta_pos $VALIDPOSFILE --valid_fasta_neg $VALIDNEGFILE \
	--prefix $PREFIX --out_dir $PROJDIR

## score validation sequences with default CNN models with default parameters
python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' \
	--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
	--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
	--train_fasta_pos $TRAINPOSFILE --train_fasta_neg $TRAINNEGFILE \
	--valid_fasta_pos $VALIDPOSFILE --valid_fasta_neg $VALIDNEGFILE \
	--prefix $PREFIX --out_dir $PROJDIR
done

## input validation files for 10-fold cross evalutation
# TESTPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_testPos.fa
# TESTNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_testNeg10x.fa
# LABEL=${CELLTYPE}_${FOLD}_training

# # score test sequences with default CNN models with default parameters
# python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' \
# --conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
# --epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
# --base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
# 	$LABEL $TRAINPOSFILE $TRAINNEGFILE $TESTPOSFILE $TESTNEGFILE 

# for MODELFILE in $PROJDIR/models/$PREFIX/*.h5; do
# FILE=$(echo $MODELFILE | sed 's/.h5/.performance.feather/g;s/\/models\//\/predictions\//g')
# if [[ ! -f "$FILE" ]]; then
# python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' \
# --train_fasta_pos $TRAINPOSFILE --train_fasta_neg $TRAINNEGFILE \
# --valid_fasta_pos $VALIDPOSFILE --valid_fasta_neg $VALIDNEGFILE \
# --out_dir $PROJDIR --model_name $MODELFILE 
# else 
# echo "Performance file exists at $(basename $FILE)."
# fi
# done 
# done



