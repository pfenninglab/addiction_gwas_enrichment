#!/bin/bash
#SBATCH -p pool1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=bam2bw
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --array=1-12
#SBATCH --error=logs/bwMean_Mouse_%A_%a.txt
#SBATCH --output=logs/bwMean_Mouse_%A_%a.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/mouse_supplement
chromSize=/home/bnphan/resources/genomes/bwa/mm10/mm10.chrom.size
cd $PROJDIR

declare -a LABEL_LIST=("" "Ctx_bulk" "Ctx_EXC" "Ctx_PV" "Ctx_SST" "Ctx_VIP" "Cpu_bulk" "Cpu_D1" "Cpu_D2" "Cpu_PV" "Cpu_SST" "NAcc_D1" "NAcc_D2")
LABEL=${LABEL_LIST[${SLURM_ARRAY_TASK_ID}]}
BW=$(ls bigWig/*${LABEL}.bw)
MEANBW=${PROJDIR}/meanBigWig/${LABEL}_mean.bw

##########################################
# create the coverage file if no bw found
if [ ! -f $MEANBW ]; 
then 
	wiggletools mean $BW | wigToBigWig stdin $chromSize $MEANBW
else 
	echo "Mean BigWig file found."; 
fi

###############################################
# compute the ATAC-seq TSS scores around genes
GTF=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/mouse_supplement/bed/gencode.vM25.annotation.gtf
BLACKLIST=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/mouse_supplement/bed/mm10.blacklist.bed
OUTFILE=${PROJDIR}/tables/${LABEL}_computeMatrix_TSS.tab.gz

if [ ! -f $OUTFILE ]; 
then 
computeMatrix reference-point \
--scoreFileName $MEANBW \
--regionsFileName $GTF \
--afterRegionStartLength 2000 \
--beforeRegionStartLength 2000 \
--metagene --referencePoint TSS \
--blackListFileName $BLACKLIST \
--samplesLabel $LABEL \
--verbose \
--outFileName $OUTFILE
else 
	echo "TSS matrix found."; 
fi




