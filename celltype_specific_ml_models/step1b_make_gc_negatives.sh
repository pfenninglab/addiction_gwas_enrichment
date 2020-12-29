#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu
#SBATCH --time 0-12
#SBATCH --job-name=genNullSeq
#SBATCH --mem=10G
#SBATCH --error=logs/genNullSeq_%A_%a_out.txt
#SBATCH --output=logs/genNullSeq_%A_%a_out.txt
#SBATCH --array=1-11

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
cd $PROJDIR

###########################################
# find the fasta file to make negatives for 
FOLD=10
FGFASTA=`ls -l $PROJDIR/FASTA/*_reproduciblePeak_summitCentered_w501.fa | awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND {print $9}'`
BGFASTA="$(sed s/_reproduciblePeak_summitCentered_w501.fa/_biasAwaySlidingGenomicBg${FOLD}x.fa/g <<<$FGFASTA | sed 's/..\/bam\///g')"

############################################################
# generate 501 bp background repository of mm10 sequences 
BGDIR=$HOME/resources/biasaway/mm10/501bp
GENOME=/home/bnphan/resources/bwa_indices/mm10/mm10.fa
if [ ! -d "$BGDIR" ]; then
	bash $HOME/resources/biasaway/create_background_repository.sh \
		-f $GENOME -s 501 -r $BGDIR
fi

# generate 501bp background sequences to BICCN huMOp cell type reproducible peaks
biasaway c --foreground $FGFASTA --nfold $FOLD \
	--deviation 2.6 --step 50 --seed 1 --winlen 100 \
	--bgdirectory $BGDIR > $BGFASTA
