options(stringsAsFactors = F, quietly = T)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

################################################
# extract the idr optimal peaks in introns and distal intergenic regions
fullardFilePath = '../gwas_enrichment/enhancers/human_brain_fansATACseq_fullard'
df = data.frame(bedFile = list.files(fullardFilePath, pattern = 'narrowPeak.gz$',recursive = T, full.names = T))
df$Region = ss(basename(df$bedFile), '_',1)
df$Celltype = ss(basename(df$bedFile), '_',2)

# keep regions enriched in Fullard files
regionKeep = c('OFC','VLPFC','DLPFC','ACC','STC','NAC','PUT')
celltypeKeep = c('NeuN+')
df = df[df$Region %in% regionKeep ,]
df = df[df$Celltype %in% celltypeKeep ,]

# make 
df$fastaFile501 = file.path('FASTA',paste0(df$Region, '_',df$Celltype,'_Fullard_IDR_optimal_peaks_501.fasta'))
df$fastaFile1000 = file.path('FASTA',paste0(df$Region, '_',df$Celltype,'_Fullard_IDR_optimal_peaks_1000.fasta'))

# read in bed peaks
peakList = GRangesList(lapply(df$bedFile, function(x){
  peaks <- sortSeqlevels( import(x))
  peaks <- sort(peaks)
  bedFile = file.path('BED',gsub('.narrowPeak.gz','.bed',basename(x)))
  export(peaks,con = bedFile)
  
  # only use standard chromosomes
  peaks = peaks[seqnames(peaks) %in% paste0('chr',c(1:22,'X','Y'))]
  peaks = setNames(peaks,paste0(seqnames(peaks),':',start(peaks),'-',end(peaks)))
  peaks = peaks[!duplicated(peaks)]
  return(peaks)
}))

###############################################
# write peak fasta sequences to file
# peakLength1 = 501 #MPRA design for 300bp oligos, 133 for barcode + adapters
peakLength1 = 501 #MPRA design for 300bp oligos, 133 for barcode + adapters
offset1 <- (peakLength1-1)/2 # make for 501 bp sequence
peakLength2 = 1000 # for CNN model predicting on 1kb sequences
offset2 <- (peakLength2-1)/2 

for(i in seq(nrow(df))){
  peaks = peakList[[i]]
  
  # only use standard chromosomes
  peaks = peaks[seqnames(peaks) %in% paste0('chr',c(1:22,'X','Y'))]
  peaks = setNames(peaks,paste0(seqnames(peaks),':',start(peaks),'-',end(peaks)))
  peaks501 = peaks1000 = peaks
  
  # for 501bp peaks
  start(peaks501) = start(peaks501) + peaks501$peak - (peakLength1 + 1)/2
  end(peaks501) =  start(peaks501) + peakLength1 -1
  
  # for 1000bp peaks
  start(peaks1000) = start(peaks1000) + peaks1000$peak - (peakLength2 + 1)/2
  end(peaks1000) =  start(peaks1000) + peakLength2 -1
  
  # get fasta sequences 
  fasta501 <- getSeq(genome, peaks501)
  fasta1000 <- getSeq(genome, peaks1000)

  # write ref/alt fasta of snps
  writeXStringSet(fasta501, df$fastaFile501[i], format = 'fasta',width = 501)
  writeXStringSet(fasta1000, df$fastaFile1000[i], format = 'fasta',width = 501)
}

save(df, peakList, file = 'rdas/fullard_idr_optimal_peaks.rda')


