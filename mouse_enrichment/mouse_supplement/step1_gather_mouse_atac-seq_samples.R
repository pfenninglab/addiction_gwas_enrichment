options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Rsubread)
library(SummarizedExperiment)
library(tidyverse)

############################
## get bulk tissue samples
DIR1 = '/projects/pfenninggroup/mouseCxStr/Aging_C57/Peaks'
folders1 = c('CTX_2_mo','CTX_3_mo','STR_2_mo','STR_3_mo')
bamFiles1 = unlist(lapply(file.path(DIR1,folders1),list.files, recursive = T, 
                  full.names = T, pattern = 'trim.merged.nodup.bam$'))
df1 = data.frame(bamFiles = bamFiles1) %>% 
  mutate(Sample = ss(basename(bamFiles), '_',1), 
         Tissue = rep(c('Ctx','Cpu'), each = 4),
         Celltype = rep(c('bulk'), each = 8), 
         Batch = 'batch1')


## Mo et al. samples
DIR2 = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/Mo2015Data'
folders2 = c('AtacExcitatoryOut','AtacPVOut','AtacVIPOut')
bamFiles2 = unlist(lapply(file.path(DIR2,folders2),list.files, recursive = T, 
                          full.names = T, pattern = '.nodup.bam$'))
df2 = data.frame(bamFiles = bamFiles2) %>% 
  mutate(Sample = ss(basename(bamFiles), '_',1), 
         Tissue = rep(c('Ctx'), each = 6),
         Celltype = rep(c('EXC','PV','VIP'), each = 2), 
         Batch = 'batch2')



## Pfenning PV-Ctx
DIR3 = '/projects/pfenninggroup/mouseCxStr/PV_INTACT_M1_CTX/Nextseq/M1_PVpos_out'
bamFiles3 = list.files(DIR3, recursive = T, full.names = T, pattern = '.nodup.bam$')
df3 = data.frame(bamFiles = bamFiles3) %>% 
  mutate(Sample = ss(basename(bamFiles), '_',1), 
         Tissue = rep(c('Ctx'), each = 3),
         Celltype = rep(c('PV'), each = 3), 
         Batch = 'batch3')


## Pfenning PV-Str
DIR4 = '/projects/pfenninggroup/mouseCxStr/parkinsons/NovaseqData/KundajePipelineOutputs/out_Saline_STR_PVpos'
bamFiles4 = list.files(DIR4, recursive = T, full.names = T, pattern = '.nodup.bam$')
df4 = data.frame(bamFiles = bamFiles4) %>% 
  mutate(Sample = ss(basename(bamFiles), '_',2), 
         Tissue = rep(c('Cpu'), each = 2),
         Celltype = rep(c('PV'), each = 2), 
         Batch = 'batch4')

## Pfenning D1,D2
DIR5 = '/projects/pfenninggroup/mouseCxStr/MotorLearning_Novaseq/ATACseqPipelineOutputs'
folders5 = list.files(DIR5, pattern = 'Naive' )
folders5 = grep('\\+',folders5, value = T)
bamFiles5 = list.files(file.path(DIR5,folders5), recursive = T, 
                       full.names = T, pattern = '.nodup.no_chrM_MT.bam$')
df5 = data.frame(bamFiles = bamFiles5) %>% 
  mutate(Sample = ss(basename(bamFiles), '_',1), 
         Tissue = rep(ss(folders5,'\\.',2), each = 2),
         Celltype = rep(ss(folders5,'\\+',1), each = 2), 
         Batch = 'batch5')

## Pfenning SST
DIR6 = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase'
folders6 = c('CortexSSTPlusYoung','StriatumSSTPlusYoung')
bamFiles6 = list.files(file.path(DIR6,folders6), recursive = T, 
                       full.names = T, pattern = '.merged.nodup.no_chrM.bam$')
bamFiles6 = bamFiles6[!grepl('glob',bamFiles6)]
df6 = data.frame(bamFiles = bamFiles6) %>% 
  mutate(Sample = ss(basename(bamFiles), '_',1), 
         Tissue = rep(c('Ctx','Cpu'), each = 2),
         Celltype = rep('SST', each = 4), 
         Batch = 'batch6')

## combine all the data
pd = do.call('rbind',list(df1, df2, df3, df4, df5, df6))
pd$Tissue = factor(pd$Tissue, levels = c('Ctx','Cpu', 'NAcc'))
pd$Celltype = factor(pd$Celltype, c('EXC','PV','SST','VIP','D1','D2','bulk'))
all(file.exists(pd$bamFiles))

saveRDS(pd, file= 'rdas/cSNAIL_phenotype_table_n31.RDS')
write.table(pd, file = 'tables/cSNAIL_phenotype_table_n31.txt', 
            sep = '\t', quote = F, row.names = F)




###################################
## make the featureCounts object
idrPeakFiles = list.files(path = '../peaks', pattern = 'narrowPeak.gz', full.names = T)
colnames = c('seqnames', 'start', 'end', 'name', 'score', 'strand',
             'signalValue', 'pValue', 'qValue', 'peak')
peakList = lapply(idrPeakFiles, read_tsv, col_names = colnames)
peakRangeList = GRangesList(lapply(peakList, GRanges))

## merge peaks and combine regions w/in 50bp
peaks = GenomicRanges::reduce( unlist(peakRangeList), min.gapwidth = 50) 
peaks2 = peaks
start(peaks2) = start(peaks2) -200
end(peaks2) = end(peaks2) + 200

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(peaks2, tssRegion=c(-5000, 5000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakRanges = as.GRanges(peakAnno)
peakRanges = setNames(peakRanges,paste0(seqnames(peakRanges),':',
                                        start(peakRanges),'-',end(peakRanges)))

gtfFile = 'tables/cSNAIL_IDR_merged-flank200_peaks.gtf'
export(as.GRanges(peakAnno), gtfFile)

## get the peak counts 
fCounts = featureCounts(files = pd$bamFile, annot.ext = gtfFile, 
                        isGTFAnnotationFile = TRUE, GTF.attrType="geneId", 
                        GTF.featureType = 'sequence_feature', 
                        useMetaFeatures = FALSE, primaryOnly = TRUE, 
                        isPairedEnd = TRUE, nthreads = 16)

colnames(fCounts$counts) = pd$Sample
rownames(fCounts$counts) = names(peakRanges)

rse = SummarizedExperiment(assays=list(counts = fCounts$counts),
                           rowRanges=peakRanges, colData=pd)
saveRDS(rse, file = 'rdas/cSNAIL_featureCounts_RangeSummarizedExperiment_n31_20200827.RDS')
