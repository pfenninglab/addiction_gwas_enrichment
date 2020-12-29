ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(SummarizedExperiment)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

DIR='/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models'
peakLength = 501 #MPRA design for 300bp oligos, 133 for barcode + adapters, 501 for ArchR compatibility

# one test set for fold cross validation
testSet =  c('chr1', 'chr2', 'chr19')

# define 8-fold cross validation chromosome hold-out scheme
validSet = list(`fold1` = c('chr6',  'chr13', 'chr21'),
                `fold2` = c('chr7',  'chr14', 'chr18'),
                `fold3` = c('chr11', 'chr17', 'chrX'),
                `fold4` = c('chr9',  'ch12'),
                `fold5` = c('chr10', 'chr8'))

#########################
# get all the cell types name
DIR='/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models'
PEAK_DIR='/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/peaks'
FASTA_SPLIT_DIR = file.path(DIR,'FASTA_CV')
BED_SPLIT_DIR = file.path(DIR,'BED_CV')

# get peak files,
bedFiles = list.files(path = PEAK_DIR, pattern = '.narrowPeak.gz', full.names = T)
bedFiles = bedFiles[!grepl('bulk',bedFiles)]
names(bedFiles) = celltypes = ss(basename(bedFiles),'.narrowPeak.gz')
idrPeakList = lapply(bedFiles, import)

df = data.frame(celltypes = rep(celltypes, each = length(validSet)), 
                fold = rep(paste0('fold', 1:length(validSet)), times = length(celltypes)))
rownames(df) = with(df, paste(celltypes, fold, sep = '_'))

df$numTestPos = NA
df$numTestNeg = NA

df$numTrainPos = NA
df$numTrainNeg = NA

df$numValidPos = NA
df$numValidNeg = NA

######################
# get all IDR peaks
for(name in celltypes){
  print(paste('Reading in fasta files of:', name))
  
  # read in the IDR fasta sequence and GC-matched negative
  posFasta = readDNAStringSet(file.path(DIR, 'FASTA', paste0(name, '_reproduciblePeak_summitCentered_w501.fa')))
  negFasta = readDNAStringSet(file.path(DIR, 'FASTA', paste0(name, '_biasAwaySlidingGenomicBg10x.fa')))
  
  # exclude sequences w/ N's
  posFasta = posFasta[!grepl('N|n',posFasta)]
  negFasta = negFasta[!grepl('N|n',negFasta)]
  
  # rename negative fasta sequence
  names(negFasta) = paste0(ss(names(negFasta), '_', 1), ':',
                           ss(names(negFasta), '_', 2), '-', 
                           as.numeric(ss(names(negFasta), '_', 2)) + 501)
  
  # exclude GC-matched negatives that are an IDR peak
  negPeak = GRanges(names(negFasta))
  start(negPeak) = start(negPeak) - 500
  end(negPeak) = end(negPeak) + 500
  indKeep = which(countOverlaps(negPeak, idrPeakList[[name]]) ==0)
  negFasta = negFasta[indKeep]
 
  # write out test sequences 
  testPosFN = file.path(FASTA_SPLIT_DIR, paste0(name, '_testPos.fa'))
  indTestPos = which(ss(names(posFasta), ':') %in% testSet)
  writeXStringSet(posFasta[indTestPos], testPosFN, format = 'fasta',width = peakLength)
  testNegFN = file.path(FASTA_SPLIT_DIR, paste0(name, '_testNeg10x.fa'))
  indTestNeg = which(ss(names(negFasta), ':') %in% testSet)
  writeXStringSet(negFasta[indTestNeg], testNegFN, format = 'fasta',width = peakLength)
  
  tmpGR = GRanges(c(names(posFasta[indTestPos]), names(negFasta[indTestNeg])))
  df1 <- data.frame(seqnames=seqnames(tmpGR), starts=start(tmpGR), ends=end(tmpGR), 
                   names = c(rep(1,length(indTestPos) ), rep(0,length(indTestNeg))))
  df1 = df1[sample(nrow(df1)),]
  write.table(file=file.path(BED_SPLIT_DIR, paste0(name, '_test.bed')), 
              df1, quote=F, sep="\t", row.names=F, col.names=F)
  
  # split test and training sequences by fold
  print(paste('Splitting test and training sequences 5-fold:', name))
  for(fold in names(validSet)){
    # # for positive sequences
    # validPosFN = file.path(FASTA_SPLIT_DIR, paste0(name, '_', fold, '_validPos.fa'))
    indValidPos = which(ss(names(posFasta), ':') %in% validSet[[fold]])
    # writeXStringSet(posFasta[indValidPos], validPosFN, format = 'fasta',width = peakLength)
    # 
    # trainPosFN = file.path(FASTA_SPLIT_DIR, paste0(name, '_', fold, '_trainPos.fa'))
    indTrainPos = which(!ss(names(posFasta), ':') %in% c(testSet,validSet[[fold]]))
    # writeXStringSet(posFasta[indTrainPos], trainPosFN, format = 'fasta',width = peakLength)
    # 
    # # for negative sequences
    # validNegFN = file.path(FASTA_SPLIT_DIR, paste0(name, '_', fold, '_validNeg10x.fa'))
    indValidNeg = which(ss(names(negFasta), ':') %in% validSet[[fold]])
    # writeXStringSet(negFasta[indValidNeg], validNegFN, format = 'fasta',width = peakLength)
    # 
    # trainNegFN = file.path(FASTA_SPLIT_DIR, paste0(name, '_', fold, '_trainNeg10x.fa'))
    indTrainNeg = which(!ss(names(negFasta), ':') %in% c(testSet,validSet[[fold]]))
    # writeXStringSet(negFasta[indTrainNeg], trainNegFN, format = 'fasta',width = peakLength)

    tmpGR2 = GRanges(c(names(posFasta[indTrainPos]), names(negFasta[indTrainNeg])))
    df2 <- data.frame(seqnames=seqnames(tmpGR2), starts=start(tmpGR2), ends=end(tmpGR2), 
                      names = c(rep(1,length(indTrainPos) ), rep(0,length(indTrainNeg) )))
    df2 = df2[sample(nrow(df2)),]
    write.table(file = file.path(BED_SPLIT_DIR, paste0(name, '_', fold, '_train.bed')), 
                df2, quote=F, sep="\t", row.names=F, col.names=F)
    
    tmpGR3 = GRanges(c(names(posFasta[indValidPos]), names(negFasta[indValidNeg])))
    df3 <- data.frame(seqnames=seqnames(tmpGR3), starts=start(tmpGR3), ends=end(tmpGR3), 
                     names = c(rep(1,length(indValidPos) ), rep(0,length(indValidNeg) )))
    df3 = df3[sample(nrow(df3)),]
    write.table(file = file.path(BED_SPLIT_DIR, paste0(name, '_', fold, '_valid.bed')), 
                df3, quote=F, sep="\t", row.names=F, col.names=F)
    
    # # make sure splitting correct
    # print(all(!duplicated(c(indTestPos, indValidPos, indTrainPos)))) 
    # print(all(!duplicated(c(indTestNeg, indValidNeg, indTrainNeg)))) 
    # 
    # # write metrics to data frame
    # df[paste(name, fold, sep = '_'), 'numTestPos'] = length(indTestPos)
    # df[paste(name, fold, sep = '_'), 'numTrainPos'] = length(indTrainPos)
    # df[paste(name, fold, sep = '_'), 'numValidPos'] = length(indValidPos)
    # df[paste(name, fold, sep = '_'), 'numTestNeg'] = length(indTestNeg)
    # df[paste(name, fold, sep = '_'), 'numTrainNeg'] = length(indTrainNeg)
    # df[paste(name, fold, sep = '_'), 'numValidNeg'] = length(indValidNeg)
  }
}

# # get the average 
# df$posToNegRate = with(df, (numTestPos/numTestNeg + 
#                             numTrainPos/numTrainNeg + 
#                             numValidPos/numValidNeg)/3)
# 
# fileName = 'tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt'
# write.table(df, file = fileName,quote = F, row.names = F, col.names = T, sep = '\t')
# 
