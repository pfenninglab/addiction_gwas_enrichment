options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(RColorBrewer)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)

####################################################
# load in annotated SNPs w/ CNN scores for each SNP
out_rds = '../gwas_enrichment/rdas/addiction_snp_cnn_scores_qval_20210413.rds'
df_snp_long = readRDS(out_rds)
df_snp_long %>% filter(tier == 'TierA') %>% pull(snp_delta_prob) %>% summary()

df_snp_long2 = df_snp_long %>% distinct(rsID, group, .keep_all = T)
snp_gr = df_snp_long %>% arrange(r2) %>% distinct(rsID, .keep_all = T) %>%
  mutate(chr = paste0('chr',chr), start = pos, end = pos, score = r2) %>% 
  column_to_rownames('rsID') %>% data.frame() %>% GRanges()

# export.bed(snp_gr, file.path('BED','addiction_gwas_all.bed'))

#lead/independent significant SNP track
snpRSE_lead = snp_gr[df_snp_long2 %>% filter(rsID == IndSigSNP) %>% pull(rsID)]
# export.bed(snpRSE_lead, file.path('BED','addiction_gwas_IndSigSNP.bed'))


#SNPs in NeuN + peaks
snpRSE_peaks = snp_gr[df_snp_long2 %>% filter(inNeuN) %>% pull(rsID)]
# export.bed(snpRSE_peaks, file.path('BED','addiction_gwas_inNeuNPeaks.bed'))

#####################################################
# export the different prioritized SNPs split by tiers
export.bed(snp_gr[df_snp_long %>% distinct(rsID, .keep_all = T) %>% 
                    filter(tier == 'TierA') %>% pull(rsID)], 
           file.path('BED','addiction_gwas_TierA_candidate.bed'))
export.bed(snp_gr[df_snp_long %>% distinct(rsID, .keep_all = T) %>% 
                    filter(tier == 'TierB') %>% pull(rsID)], 
           file.path('BED','addiction_gwas_TierB_candidate.bed'))
export.bed(snp_gr[df_snp_long %>% distinct(rsID, .keep_all = T) %>% 
                    filter(tier == 'TierC') %>% pull(rsID)], 
           file.path('BED','addiction_gwas_TierC_candidate.bed'))

########################################
# export the Tier A SNP split by group 
tmp = df_snp_long %>% filter(tier == 'TierA') %>% distinct(rsID, group, .keep_all = T) %>% 
  mutate(chr = paste0('chr',chr), start = pos, name = rsID, end = pos, score = snp_delta_prob)
tmpList = split(tmp, tmp$group)
sapply(tmpList, nrow)

lapply(names(tmpList), function(name){
  export.bed(GRanges(tmpList[[name]]), 
             file.path('BED',paste0('addiction_gwas_CNN_scored_',name,'.bed')))
})

##########################################
# get loci to plot around candidate SNPs
snpBlock =  df_snp_long %>% arrange(chr, pos) %>% filter(tier == 'TierA') %>% 
  distinct(GenomicLocus) %>% mutate(name = GenomicLocus, value = GenomicLocus) %>% 
  select(value, name) %>% deframe() %>% GRanges()
start(snpBlock) = floor(start(snpBlock)/10^3 - 1)*10^3+1 # round to nearest 1kb
end(snpBlock) = ceiling(end(snpBlock)/10^3 + 1)*10^3 # round to nearest 1kb
snpBlock = GenomicRanges::reduce(snpBlock)
export.bed(snpBlock, file.path('BED','snp_candidate_snpBlock.bed'))


###########################################################
# filter GTF file just around snpBlock, for fast plotting
gtf = import('gene_tracks/gencode.v32.basic.annotation.gtf.gz')

# filter out transcript_type and transcript_support_level
# https://www.gencodegenes.org/pages/data_format.html
gtf = gtf[gtf$transcript_support_level %in% c(1, 2)]
gtf = gtf[gtf$transcript_type %in% c('protein_coding', 'retained_intron', 
                                     'lncRNA','polymorphic_pseudogene')]

# round to nearest 10kb
snpBlock2 = snpBlock
start(snpBlock2) = floor(start(snpBlock2)/5/10^4 - 1)*5*10^4+1 - 10^5 
# round to nearest 10kb
end(snpBlock2) = ceiling(end(snpBlock2)/5/10^4 + 1)*5*10^4 + 10^5
oo = findOverlaps(query = gtf, subject = snpBlock2)
gtf_filtered = gtf[unique(queryHits(oo))]
export(gtf_filtered, 'gene_tracks/gencode.v32.basic.annotation.filtered.gtf')



