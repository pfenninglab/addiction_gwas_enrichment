#!/bin/bash
#SBATCH -p pool1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=bam2bw
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --array=1-31
#SBATCH --error=logs/bam2bw_Mouse_%A_%a.txt
#SBATCH --output=logs/bam2bw_Mouse_%A_%a.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/mouse_supplement

cd $PROJDIR
((ID=$SLURM_ARRAY_TASK_ID + 1))
BAMFILE=$(awk -v j="$ID" 'NR==j {print $1}' $PROJDIR/tables/cSNAIL_phenotype_table_n31.txt | tr -d '[:space:]')
SAMPLE=$(awk -v j="$ID" 'NR==j {print $2}' $PROJDIR/tables/cSNAIL_phenotype_table_n31.txt | tr -d '[:space:]')
TISSUE=$(awk -v j="$ID" 'NR==j {print $3}' $PROJDIR/tables/cSNAIL_phenotype_table_n31.txt | tr -d '[:space:]')
CELLTYPE=$(awk -v j="$ID" 'NR==j {print $4}' $PROJDIR/tables/cSNAIL_phenotype_table_n31.txt | tr -d '[:space:]')
BWFILE=${PROJDIR}/bigWig/${SAMPLE}_${TISSUE}_${CELLTYPE}.bw


###############################
# index bamfile if no bai file
if [ ! -f $(echo ${BAMFILE} | sed 's/.bam/.bam.bai/g') ]; 
then 
	samtools index $BAMFILE; 
else 
	echo "BAM index found."; 
fi

##########################################
# create the coverage file if no bw found
if [ ! -f $BWFILE ]; 
then 
	bamCoverage --bam $BAMFILE --outFileName $BWFILE --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX chrY --extendReads
else 
	echo "bigWig file found."; 
fi




