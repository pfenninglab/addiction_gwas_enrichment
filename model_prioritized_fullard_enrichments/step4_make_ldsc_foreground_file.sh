## for CNN Scored backgrounds ## 
OUTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/model_prioritized_fullard_enrichments
BG=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/HoneyBadger/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
MERGED_BG=BED/mergedHoneyBadger_Fullard_CNN_scored.bed.gz
BG_NAME='HoneyBadger_Fullard_CNN_scored'
ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Fullard_CNNscored_EUR2
cd $OUTDIR
mkdir -p $ANNOTDIR
# rm $ANNOTDIR/*
# merge peaks with honeybadger2 DHS background
# find 'BED' -name '*cnn_scored_Fullard_peaks.bed' -exec cat {} \; > BED/tmp.bed
# zcat BED/*cnnV3_scored_Fullard_peaks.bed.gz $BG > BED/tmp.bed
# cut -f 1-3 BED/tmp.bed > BED/merged_honeybadger.tmp.bed
# sort -k1,1 -k2,2n BED/merged_honeybadger.tmp.bed > BED/merged_honeybadger.sorted.tmp.bed 
# bedtools merge -i BED/merged_honeybadger.sorted.tmp.bed | gzip > $MERGED_BG
# rm BED/tmp.bed BED/*.tmp.bed

# # annotated the backgrounds
# sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh \
# 	-i $MERGED_BG -n $BG_NAME -o $ANNOTDIR

CTS_FN=$OUTDIR/Fullard_cnn_scored_peaks.ldcts; > $CTS_FN
# annotate the foregrounds
for j in {1..36};
do
	NAME=$(awk -v j="$j" 'NR==j {print $1}' Fullard_cnn_scored_peaks.txt | tr -d '[:space:]')
	INPUT=$(awk -v j="$j" 'NR==j {print $2}' Fullard_cnn_scored_peaks.txt | tr -d '[:space:]')
	# sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i ${INPUT} -n ${NAME} -o $ANNOTDIR
	echo -e "${NAME}\t${ANNOTDIR}/${NAME}.,${ANNOTDIR}/${BG_NAME}." >> $CTS_FN
done

# ################################
# ## for SVM scored foregrounds ##
# OUTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/model_prioritized_fullard_enrichments/ldsc_enrichments/Fullard_svm_scored_enrichments
# BG=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/HoneyBadger/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed
# MERGED_BG=BED/mergedHoneyBadger_Fullard_SVM_scored.bed
# BG_NAME='HoneyBadger_Fullard_SVM_scored'
# ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Fullard_SVMscored_EUR

# # merge peaks with honeybadger2 DHS background
# find 'BED' -name '*svm_scored_Fullard_peaks.bed' -exec cat {} \; > BED/tmp.bed
# cat  BED/tmp.bed $BG > BED/tmp2.bed
# cut -f 1-3 BED/tmp2.bed > BED/merged_honeybadger.tmp.bed
# sort -k1,1 -k2,2n BED/merged_honeybadger.tmp.bed > BED/merged_honeybadger.sorted.tmp.bed 
# bedtools merge -i BED/merged_honeybadger.sorted.tmp.bed  > $MERGED_BG
# rm BED/tmp.bed BED/tmp2.bed BED/*.tmp.bed

# # annotated the backgrounds
# sbatch --partition pfen3 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh \
# 	-i $MERGED_BG \
# 	-n $BG_NAME \
# 	-o $ANNOTDIR

# # annotate the foregrounds
# for j in {1..19};
# do
# 	NAME=$(awk -v j="$j" 'NR==j {print $1}' Fullard_svm_scored_peaks.txt | tr -d '[:space:]')
# 	INPUT=$(awk -v j="$j" 'NR==j {print $2}' Fullard_svm_scored_peaks.txt | tr -d '[:space:]')
# 	sbatch --partition pfen3 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i ${INPUT} -n ${NAME} -o $ANNOTDIR
# done



# #################################################
# ## for combined SVM and CNN scored foregrounds ##
# OUTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/model_prioritized_fullard_enrichments/ldsc_enrichments/Fullard_both_scored_enrichments
# BG=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/HoneyBadger/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed
# MERGED_BG=BED/mergedHoneyBadger_Fullard_both_scored.bed
# BG_NAME='HoneyBadger_Fullard_both_scored'
# ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/Fullard_BOTHscored_EUR

# # merge peaks with honeybadger2 DHS background
# find 'BED' -name '*both_scored_Fullard_peaks.bed' -exec cat {} \; > BED/tmp.bed
# cat  BED/tmp.bed $BG > BED/tmp2.bed
# cut -f 1-3 BED/tmp2.bed > BED/merged_honeybadger.tmp.bed
# sort -k1,1 -k2,2n BED/merged_honeybadger.tmp.bed > BED/merged_honeybadger.sorted.tmp.bed 
# bedtools merge -i BED/merged_honeybadger.sorted.tmp.bed  > $MERGED_BG
# rm BED/tmp.bed BED/tmp2.bed BED/*.tmp.bed

# # annotated the backgrounds
# sbatch --partition pfen3 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh \
# 	-i $MERGED_BG \
# 	-n $BG_NAME \
# 	-o $ANNOTDIR

# # annotate the foregrounds
# for j in {1..19};
# do
# 	NAME=$(awk -v j="$j" 'NR==j {print $1}' Fullard_both_scored_peaks.txt | tr -d '[:space:]')
# 	INPUT=$(awk -v j="$j" 'NR==j {print $2}' Fullard_both_scored_peaks.txt | tr -d '[:space:]')
# 	sbatch --partition pfen3 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i $INPUT -n $NAME -o $ANNOTDIR
# done




