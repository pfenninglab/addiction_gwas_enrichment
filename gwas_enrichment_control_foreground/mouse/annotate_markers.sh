#!/bin/bash

#ANNOTATE CELL TYPE FOREGROUNDS
for file in *3.bed;
do
  sbatch --partition pfen_bigmem /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i $file -n "${file::-12}" -o /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Controls_EUR/mouse
done

#ANNOTATE BACKGROUNDS

#HONEYBADGER2 + ALL FOREGROUNDS MERGED
sbatch --partition pfen_bigmem /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/foreground/controls/mouse/mouseDHS_merged.bed -n MouseDHS -o /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Controls_EUR/mouse
