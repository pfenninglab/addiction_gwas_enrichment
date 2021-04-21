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
group = c('Ctx-EXC','Ctx-PV','Ctx-SST','Ctx-VIP',
          'Str-PV','Str-SST','Cpu-D1','Cpu-D2', 'Nac-D1','Nac-D2')
models = read_tsv('../celltype_specific_ml_models/tables/SingleTask_IDR_bestModels_CNN_performance.tsv')
models$sample2 = ss(models$sample, '_fold')
# models$celltypes = factor(ss(ss(models$sample,'_',2), 'pos'), celltypes)
# models$tissue = factor(ss(models$sample,'_',3), tissue)
# models$tissue2 = ifelse(models$tissue== 'Ctx','Ctx','Str')
models$group = factor(paste(models$tissue, models$celltype, sep = '-'))

# get model names and folds, collapse folds/samples across tissue/celltype
splitInd = split(seq(nrow(models)), models$group)
matchList = c(`ACC_NeuN+` = 'Ctx-EXC',
                 `DLPFC_NeuN+` = 'Ctx-EXC',
                 `NAC_NeuN+` = 'Nac-D1',
                 `NAC_NeuN+` = 'Nac-D2',
                 `OFC_NeuN+` = 'Ctx-EXC',
                 `PUT_NeuN+` = 'Cpu-D1',
                 `PUT_NeuN+` = 'Cpu-D2',
                 `STC_NeuN+` = 'Ctx-EXC',
                 `VLPFC_NeuN+` = 'Ctx-EXC')
all(matchList %in% names(splitInd))
cutOffList =  1-2^(-seq(1:4))
names(cutOffList) = cutOffList

######################################
# read in CNN enhancer model scores  #
for(ind in seq_along(matchList)){
  label = paste0(names(matchList),'_Fullard_IDR_optimal_peaks')[ind]
  print(paste('Reading in',label, '.'))
  df1 = data.table::fread(file.path('tables','cnn_predictions',paste0(label, '_v3.tsv.gz'))) %>%
    left_join(models, by = 'model') %>% 
    pivot_wider(names_from = group, id_cols = 'label', values_from = 'y_pred_score', 
                values_fn = mean) %>% column_to_rownames('label')
  gr = GRanges(rownames(df1))
  df2 = df1[matchList[ind]]
  
    peakList = lapply(cutOffList,  function(cutOff) gr[df2[,1] >= cutOff])
    names(peakList) = paste0(matchList[ind], '_cutoff',cutOffList)
    print(lengths(peakList))
    
    for(name in names(peakList)){
      print(paste('processing:', name))
      bedFile = file.path('BED',paste(gsub('NeuN\\+_Fullard_IDR_optimal_peaks', name, label),
                                      'cnnV3_scored_Fullard_peaks.bed.gz', sep = '_'))
      export(peakList[[name]], con = bedFile)
  }
}
bedFiles = list.files('BED',pattern = 'cnnV3_scored_Fullard_peaks.bed.gz', full.names = T)
label = ss(basename(bedFiles), '_cnnV3_scored_Fullard_peaks.bed.gz')
df3 = data.frame(label = label,
                 bedFiles = bedFiles)
write.table(df3, file = 'Fullard_cnn_scored_peaks.txt',quote = F, sep = '\t',
            row.names = F, col.names = F)




