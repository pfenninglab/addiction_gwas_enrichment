###############################
# Step 1. make the track files
# make_tracks_file --trackFiles $OFC $VLPFC $DLPFC $STC $PUT $NAC $ALL $LEAD $CAUSAL $PEAK $CANDIDATE $GTF -o $TRACKS

DIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/mouse_supplement
mkdir -p $DIR/candidate_plots
cd $DIR/candidate_plots

# SNP blocks to plot
BLOCK1=$DIR/bed/cSNAIL_markerGene_regions.bed
TRACKS20=$DIR/config/config_meanBW_cSNAIL20.ini
TRACKS40=$DIR/config/config_meanBW_cSNAIL40.ini
TRACKS80=$DIR/config/config_meanBW_cSNAIL80.ini

# tall track for figure
pyGenomeTracks --tracks $TRACKS20 --BED $BLOCK1 --fontSize 4 --width 4 --dpi 1000 --trackLabelFraction 0.1 --outFileName tall_meanBW_cSNAIL_marker20.pdf

# tall track for figure
pyGenomeTracks --tracks $TRACKS40 --BED $BLOCK1 --fontSize 4 --width 4 --dpi 1000 --trackLabelFraction 0.1 --outFileName tall_meanBW_cSNAIL_marker40.pdf

# tall track for figure
pyGenomeTracks --tracks $TRACKS80 --BED $BLOCK1 --fontSize 4 --width 4 --dpi 1000 --trackLabelFraction 0.1 --outFileName tall_meanBW_cSNAIL_marker80.pdf
