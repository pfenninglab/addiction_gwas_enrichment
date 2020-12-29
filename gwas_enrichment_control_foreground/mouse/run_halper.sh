#!/bin/bash

declare -a arr=("kidney_female_atac_" "liver_female_atac_" "lung_male_atac_")

for i in "${arr[@]}"
do
    python ~/halLiftover-postprocessing/orthologFind.py -max_len 1000 -min_len 50 -max_frac 1.0 -min_frac 0.1 -protect_dist 5 -qFile "${i}labeled_mm10.narrowPeak" -tFile "${i}labeled_halLiftover_hg38.bed" -sFile "${i}summits_halLiftover_hg38.bed" -oFile "${i}halper.bed" -mult_keepone
done

for file in *halper.bed;
do
  sort -k 1,1 -k 2,2n $file | awk '{print $1, $2, $3}' | perl -p -i -e 's/ /\t/g' > "${file::-4}3.bed"
done
