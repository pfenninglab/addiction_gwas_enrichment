
###############################
# 14 cm x 11cm
cd /projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/model_prioritized_fullard_enrichments/genomeTracks
# SNP blocks to plot
BLOCK=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/model_prioritized_fullard_enrichments/BED/plot_region_for_pyGenomeTracks.bed
TRACKS=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/model_prioritized_fullard_enrichments/config/cortical_pyGenomeTrack_BICCN_huMOp.ini

# tall track for figure
pyGenomeTracks --tracks $TRACKS --BED $BLOCK --fontSize 7 --width 13.97 --dpi 1000 --trackLabelFraction 0.12 --outFileName model_prioritized_fullard_enrichments_BICCN_huMOp_SNARE-Seq2_EXC-PV-VIP.pdf











