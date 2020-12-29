options(stringsAsFactors = F, quietly = T)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(GenomicRanges)
library(rtracklayer)
library(feather)
library(tidyverse)
library(pracma)
load(file = 'rdas/fullard_idr_optimal_peaks.rda')
df$name = names(peakList) = with(df, paste(Region,Celltype, sep = '_'))
df = df[which(df$Celltype =='NeuN+'), ]

# read in the CNN models
celltypes = c('EXC','PV','SST','VIP','D1','D2')
tissue = c('Ctx','Cpu','Nac')
group = c('Ctx-EXC','Ctx-PV','Ctx-SST','Ctx-VIP','Str-PV','Str-SST','Str-D1','Str-D2')
models = read_feather('../celltype_specific_ml_models/rdas/SingleTask_IDR_peakPrediction_CNN_performance.feather')
models$sample2 = ss(models$sample, '_fold')
models$celltypes = factor(ss(ss(models$sample,'_',2), 'pos'), celltypes)
models$tissue = factor(ss(models$sample,'_',3), tissue)
models$tissue2 = ifelse(models$tissue== 'Ctx','Ctx','Str')
models$group = factor(paste(models$tissue2, models$celltypes, sep = '-'), group)

# get model names and folds, collapse folds/samples across tissue/celltype
splitInd = split(seq(nrow(models)), models$group)

StriatumSamples = c('PUT','NAC')

cutOffList =  1-2^(-seq(1:6))

######################################
# read in CNN enhancer model scores  #
for(label in gsub('_501.fasta', '',basename(df$fastaFile501))){
  print(paste('Reading in',label, '.'))
  df1 = read.delim(file.path('tables','cnn_predictions',paste0(label, '.tsv.gz')),
                   row.names = 1)[-c(1:2),] %>% apply(1, as.numeric)
  df2 = do.call('cbind',lapply(splitInd, function(i) colMeans(df1[i,]))) %>% data.frame()
  gr = GRanges(rownames(df2))
  
  if(grepl(paste(StriatumSamples,collapse = '|'), label)) {
    whichCols = 5:8
  } else {
    whichCols = 1:4
  }
    for(cutOff in cutOffList){
      print(paste('Using cutoff:', cutOff))

      peakList = apply(df2, 2, function(x) gr[x >= cutOff])[whichCols]
      print(lengths(peakList))
      for(name in names(peakList)){
        bedFile = file.path('BED',paste(gsub('NeuN\\+_Fullard_IDR_optimal_peaks', name, label), 
                                        paste0('cutOff',cutOff), 'cnn_scored_Fullard_peaks.bed', sep = '_'))
        export(peakList[[name]],con = bedFile)
    }
  }
}
bedFiles = list.files('BED',pattern = 'cnn_scored_Fullard_peaks.bed', full.names = T)
label = paste(ss(basename(bedFiles), '_'), 
              ss(ss(basename(bedFiles), '_',2), '\\.',2),
              ss(ss(basename(bedFiles), '_',3), 'cutOff',2), sep = '_')
df3 = data.frame(label = label,
                 bedFiles = bedFiles)
write.table(df3, file = 'Fullard_cnn_scored_peaks.txt',quote = F, sep = '\t',
            row.names = F, col.names = F)




