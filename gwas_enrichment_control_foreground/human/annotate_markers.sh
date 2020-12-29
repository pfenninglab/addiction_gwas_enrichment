#!/bin/bash


#ANNOTATE CELL TYPE FOREGROUNDS
for file in *.bed;
do
sbatch --partition pfen_bigmem /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i ${file} -n "${file::-4}" -o /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Controls_EUR/human
#done

#ANNOTATE BACKGROUNDS

#HONEYBADGER2 + FOREGROUNDS MERGED
sbatch --partition pfen_bigmem /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/foreground/controls/human/honeybadger2_merged.bed -n Honeybadger -o /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Controls_EUR/human
