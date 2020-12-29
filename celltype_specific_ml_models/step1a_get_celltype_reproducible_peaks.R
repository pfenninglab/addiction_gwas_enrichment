options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(umap)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10.masked)
genome <- BSgenome.Mmusculus.UCSC.mm10.masked

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakLength = 501 #MPRA design for 300bp oligos, 133 for barcode + adapters, 501 for ArchR compatibility

DIR='/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models'
PEAK_DIR='/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/peaks'

# get peak files,
bedFiles = list.files(path = PEAK_DIR, pattern = 'narrowPeak.gz', full.names = T)
names(bedFiles) = ss(basename(bedFiles),'.narrowPeak.gz')
bedFiles = bedFiles[!grepl('bulk',bedFiles)]

# annotate peaks, and select intronic/distal intergenic
summitList = GRangesList(lapply(bedFiles, function(x){
  peaks = import(x)
  # get location of idr peak summits
  start(peaks) = start(peaks) + peaks$peak - (peakLength + 1)/2
  end(peaks) = start(peaks) + peakLength -1
  peaks$peak = start(peaks) + (peakLength + 1)/2
  peakRanges = as.GRanges(annotatePeak(peaks, tssRegion=c(-5000, 5000), TxDb=txdb, annoDb="org.Mm.eg.db"))
  peakRanges = setNames(peakRanges,paste0(seqnames(peakRanges),':',start(peakRanges),'-',end(peakRanges)))
  # filter out exons, promoters
  filteredRanges = peakRanges[ss(peakRanges$annotation, '\\(',1) %in% c('Distal Intergenic','Intron')] #distal intergenic
  return(filteredRanges)
}))


system('mkdir -p FASTA')
for(lab in names(summitList)){
  tmpGR = summitList[[lab]]
  offset <- 20 # make for 41 bp sequence
  
  # get fasta sequences for effect allele
  fasta <- getSeq(genome,seqnames(tmpGR),start(tmpGR),end(tmpGR))
  names(fasta) = paste0(seqnames(tmpGR),':',start(tmpGR),'-',end(tmpGR))
  
  # write fasta sequence
  fastaFile = paste0('FASTA/',lab,'_reproduciblePeak_summitCentered_w',peakLength,'.fa')
  writeXStringSet(fasta, fastaFile, format = 'fasta',width = peakLength)
}





