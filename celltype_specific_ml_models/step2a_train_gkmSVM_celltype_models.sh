#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen_bigmem
##SBATCH --time=1-0
#SBATCH --job-name=svm
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=25G
#SBATCH --array=31-55
#SBATCH --error=logs/lsgkm_%A_%a_output.txt
#SBATCH --output=logs/lsgkm_%A_%a_output.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
cd $PROJDIR

# cell type and fold index
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $1}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`
FOLD=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $2}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`

# input training files 
TRAINPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainPos.fa
TRAINNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainNeg10x.fa

# input validation files for k-fold cross evalutation
VALIDPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validPos.fa
VALIDNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validNeg10x.fa

LABEL=${CELLTYPE}_${FOLD}

# parameters to tune SVM
CPARAM=1 # regularization param
WPARAM=10 # weight of the positives

# train SVM model
mkdir -p $PROJDIR/svm/${LABEL}/
OUTPREFIX=$PROJDIR/svm/${LABEL}/${LABEL}_c${CPARAM}_w${CPARAM}
if [ ! -f "${OUTPREFIX}.model.txt" ]; then
	gkmtrain -t 4 -l 11 -k 7 -d 3 -c ${CPARAM} -M 50 -H 50 -w ${WPARAM} \
	-m 10000 -s -T 4 -e 0.001 $TRAINPOSFILE $TRAINNEGFILE $OUTPREFIX
fi

# pedict test set sequences
mkdir -p $PROJDIR/predictions/svm/${LABEL}/
PREDICTBASE=$PROJDIR/predictions/svm/${LABEL}/$LABEL_c${CPARAM}_w${CPARAM}
if [ ! -f "${PREDICTBASE}_gkmpredict_validPos.txt" ]; then
	gkmpredict -T 4 $VALIDPOSFILE ${OUTPREFIX}.model.txt ${PREDICTBASE}_gkmpredict_validPos.txt
	gkmpredict -T 4 $VALIDNEGFILE ${OUTPREFIX}.model.txt ${PREDICTBASE}_gkmpredict_validNeg.txt
fi
