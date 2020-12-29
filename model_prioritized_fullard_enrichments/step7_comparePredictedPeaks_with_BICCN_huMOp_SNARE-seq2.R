ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(GenometriCorr)
library(parallel)

enrichFisher = function(peak_pred, markerPeaks, universe){
  fisherMat = do.call('rbind',lapply(markerPeaks, function(markers){
    test = fisher.test(table(countOverlaps(universe, peak_pred) > 0, 
                       countOverlaps(universe, markers) > 0))
    return(data.frame(p.value = test$p.value,
                      OR = test$estimate,
                      OR.conf.int.max = max(test$conf.int),
                      OR.conf.int.min = min(test$conf.int)))
  }))
  return(fisherMat)
}

peakPrecision = function(peak_true, peak_pred){
  TP = sum(countOverlaps(peak_pred, peak_true)>=1)
  FP = sum(countOverlaps(peak_pred, peak_true)==0)
  precision = TP/(TP+ FP)
}

##################################
# get SNARE-seq2 reproducible peaks
ArchRDir ="/projects/pfenninggroup/singleCell/BICCN_human_M1_SNARE-Seq2/ArchR_analyses/ArchR_BICCN_huMOp_SNARE-Seq2_level1"
SNAREseqTracks = list(
  EXC = file.path(ArchRDir,'PeakCalls','Glutamatergic-reproduciblePeaks.gr.rds'),
  PV = file.path(ArchRDir,'PeakCalls','PVALB-reproduciblePeaks.gr.rds'),
  VIP = file.path(ArchRDir,'PeakCalls','VIP-reproduciblePeaks.gr.rds'))

peakList = GRangesList(lapply(SNAREseqTracks, readRDS))
lapply(names(SNAREseqTracks), function(x){
  gr = sort(peakList[[x]])
  fileName = file.path('BED', paste0('BICCN_human_M1_SNARE-Seq2_',x,'.bed'))
  export(gr, fileName)
  })

##################################
# get Lake et al. cell type specific peaks
lakeDIR ="/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/gwas_enrichment/enhancers"
lakeTracks = list(
  Ex = file.path(lakeDIR,'human_occ_scTHS_lake_main','GSE97887_occ1_MAINSPLITS.Ex.diffPeaks.bed'),
  In = file.path(lakeDIR,'human_occ_scTHS_lake_main','GSE97887_occ1_MAINSPLITS.In.diffPeaks.bed'))
lakePeakList = GRangesList(lapply(lakeTracks, import))

######################################################
## get model-predicted peaks from NeuN+ open chromatin
predictedTracks = list.files(path = 'BED', full.names = TRUE, 
                             pattern = '_both_scored_Fullard_peaks.bed')
predictedTracks = predictedTracks[grepl('DLPFC|OFC|STC|VLPFC|ACC', predictedTracks)]
predictedPeaks = GRangesList(lapply(predictedTracks, import))
names(predictedPeaks) = gsub('_both_scored_Fullard_peaks.bed', '',basename(predictedTracks))

#####################################################
## calculate overlap enrichment of ML-predicted peaks
universe = GenomicRanges::reduce(unlist(c(predictedPeaks, peakList)))
enrichMat = do.call('rbind',lapply(predictedPeaks, enrichFisher, 
                                   markerPeaks = peakList,
                                   universe = universe))
enrichMat$predictedGroup = ss(rownames(enrichMat),'\\.', 1)
enrichMat$celltypeGroup = ss(rownames(enrichMat),'\\.', 2)

# filter for positive enrichment and FDR < 0.05
enrichMat %>% mutate(FDR = p.adjust(p.value, 'BH')) %>% filter(FDR < 0.05, OR> 1)
mapply(peakPrecision, peak_true =peakList, peak_pred = mlPeakList)


## calculate overlap w/ GenometriCorr
permut.number <- 0

enrichMat = do.call('rbind',lapply(peakList,function(peak){ 
  outList =parallel::mclapply(predictedPeaks, GenometriCorrelation,reference = peak, awhole.only = TRUE,
                         chromosomes.to.include.in.awhole = paste0('chr',c(1:22,'X')),
                         permut.number = permut.number, keep.distributions = TRUE, 
                         showProgressBar = FALSE, mc.cores = parallel::detectCores()/2)
  out = do.call('rbind', outList)
  return(out)
  }))
    
GenometriCorrelation( predictedPeaks[[1]], peakList[[1]], awhole.only = TRUE,
                 chromosomes.to.include.in.awhole = paste0('chr',c(1:22,'X')),
                 permut.number = permut.number, keep.distributions = TRUE, 
                 showProgressBar = FALSE)


############################################
## prepare tracks for pyGenomeTrack plots


# filter GTF file just around plotRegion, for fast plotting
gtf = import('../genomeTracks/gene_tracks/gencode.v32.basic.annotation.gtf.gz')

# filter out transcript_type and transcript_support_level
# https://www.gencodegenes.org/pages/data_format.html
gtf = gtf[gtf$transcript_support_level %in% c(1, 2)]
gtf = gtf[gtf$transcript_type %in% c('protein_coding', 'retained_intron', 
                                     'lncRNA','polymorphic_pseudogene')]

###########################################
# get loci to plot around several gene loci
genes = c( 'RBFOX3', 'CAMK2A','BDNF','TCF4', 'SLC17A7', 'GAD1', 'PVALB', 'SST')
plotRegion = reduce(gtf[gtf$gene_name %in% genes,])
start(plotRegion) = floor(start(plotRegion)/10^5 - 1)*10^5+1 # round to nearest 10kb
end(plotRegion) = ceiling(end(plotRegion)/10^5 + 1)*10^5 # round to nearest 10kb
plotRegion$score = 0
export.bed(plotRegion, file.path('BED','plot_region_for_pyGenomeTracks.bed'))

start(plotRegion) = start(plotRegion) - 5* 10^5
end(plotRegion) = end(plotRegion) +  5* 10^5
oo = findOverlaps(query = gtf, subject = plotRegion)
gtf_filtered = gtf[unique(queryHits(oo))]
export(gtf_filtered, 'gene_tracks/gencode.v32.basic.annotation.filtered.gtf')





