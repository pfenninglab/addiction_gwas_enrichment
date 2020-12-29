cd /projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/gwas_enrichment/pyGenomeTracks

######################
# plot the LD block 
pyGenomeTracks --tracks Fullard_bnp.ini \
	--BED ../tables/fuma_addiction_loci.bed \
	--outFileName loci_plots/addiction_Fullard_peaks_ldBlock.pdf \
	--fontSize 10 --height 24 --width 32 --trackLabelFraction .1 # 15cm x 20cm ~ 6in x 8in

# mouse tracks loci plots
pyGenomeTracks --tracks Mouse_tracks.ini \
	--BED ../tables/fuma_addiction_loci.bed \
	--outFileName loci_plots/addiction_Mouse_peaks_ldBlock.pdf \
	--fontSize 10 --height 15 --width 32 --trackLabelFraction .1 # 15cm x 20cm ~ 6in x 8in

#############################################
# plot the 20kb region around interesting snps
pyGenomeTracks --tracks Fullard_bnp.ini \
	--BED ../tables/fuma_addiction_snps_interesting.bed \
	--outFileName snp_plots/addiction_Fullard_peaks_snp.pdf \
	--fontSize 10 --height 15 --width 8 --trackLabelFraction .1 # 15cm x 20cm ~ 6in x 8in


# mouse tracks snp plots
pyGenomeTracks --tracks Mouse_tracks.ini \
	--BED ../tables/fuma_addiction_snps_interesting.bed \
	--outFileName snp_plots/addiction_Mouse_peaks_snps.pdf \
	--fontSize 10 --height 15 --width 32 --trackLabelFraction .1 # 15cm x 20cm ~ 6in x 8in







#############################################
# plot the 20kb region around SUFU snp rs10786689, chr10:102596602
pyGenomeTracks --tracks Fullard_box.ini \
	--region chr10:102595602-102597602 \
	--outFileName snp_plots/addiction_Fullard_peaks_rs10786689.pdf \
	--fontSize 8 --height 15 --width 12 --trackLabelFraction .2 # 15cm x 20cm ~ 6in x 8in


# mouse tracks snp plots
pyGenomeTracks --tracks Mouse_tracks.ini \
	--region chr10:102595602-102597602 \
	--outFileName snp_plots/addiction_Mouse_peaks_rs10786689.pdf \
	--fontSize 8 --height 10 --width 12 --trackLabelFraction .3 # 15cm x 20cm ~ 6in x 8in






