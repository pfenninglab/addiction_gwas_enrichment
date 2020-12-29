###############################
# Step 1. make the track files
# make_tracks_file --trackFiles $OFC $VLPFC $DLPFC $STC $PUT $NAC $ALL $LEAD $CAUSAL $PEAK $CANDIDATE $GTF -o $TRACKS

DIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/genomeTracks
cd $DIR/candidate_plots

# SNP blocks to plot
BLOCK1=$DIR/BED/ctx_candidate_snpBlock.bed
TRACKS1=$DIR/config/candidate_snps_track_gtf_NeuN+_mmOrth_Ctx.ini

# tall track for figure
pyGenomeTracks --tracks $TRACKS1 --BED $BLOCK1 --fontSize 7 --width 14 --dpi 1000 --trackLabelFraction 0.12 --outFileName addiction_snp_ctx_candidate.pdf

# wide track for ppt
pyGenomeTracks --tracks $TRACKS1 --BED $BLOCK1 --fontSize 9 --width 23 --dpi 1000 --trackLabelFraction 0.12 --outFileName addiction_snp_ctx_candidate_wide.pdf


###############################
# SNP blocks to plot STR snps
BLOCK2=$DIR/BED/str_candidate_snpBlock.bed
TRACKS2=$DIR/config/candidate_snps_track_gtf_NeuN+_mmOrth_Str.ini


# tall track for figure
pyGenomeTracks --tracks $TRACKS2 --BED $BLOCK2 --fontSize 7 --width 14 --dpi 1000 --trackLabelFraction 0.12 --outFileName addiction_snp_str_candidate.pdf

# wide track for ppt
pyGenomeTracks --tracks $TRACKS2 --BED $BLOCK2 --fontSize 9 --width 23 --dpi 1000 --trackLabelFraction 0.12 --outFileName addiction_snp_str_candidate_wide.pdf


###############################
# SNP blocks to plot 
BLOCK3=$DIR/BED/snp_candidate_snpBlock.bed
TRACKS3=$DIR/config/candidate_snps_track_gtf_NeuN+_mmOrth.ini


# tall track for figure
pyGenomeTracks --tracks $TRACKS3 --BED $BLOCK3 --fontSize 7 --width 14 --dpi 1000 --trackLabelFraction 0.12 --outFileName addiction_snp_candidate.pdf

# wide track for ppt
pyGenomeTracks --tracks $TRACKS3 --BED $BLOCK3 --fontSize 9 --width 23 --dpi 1000 --trackLabelFraction 0.12 --outFileName addiction_snp_candidate_wide.pdf









