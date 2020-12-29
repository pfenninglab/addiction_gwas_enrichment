options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(RColorBrewer)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)

######################
# load in annotated SNPs
load('../gwas_enrichment/rdas/addiction_snp_withCausalDB_Fullard_overlap_rse.rda')
rowData(snpRSE) = snpDat2 = snpDat

# load in the SNP zscores
load('../celltype_specific_ml_models/rdas/addiction_snp_cnn_score_20200727.rda')
celltypes = c('Ctx-EXC','Str-D1', 'Str-D2')
rowData(snpRSE)$effect = df_effect[rownames(rowData(snpRSE)),celltypes]
df_effect = df_effect[,celltypes]

outFile='tables/addiction_snp_cnn_effectAllele_fullardPeaks_causalDB.xlsx'
# writexl::write_xlsx(cbind(snpDat, df_effect[rownames(rowData(snpRSE)),]),  outFile)
outFile2='tables/addiction_snp_cnn_nonEffAllele_fullardPeaks_causalDB.xlsx'
# writexl::write_xlsx(cbind(snpDat, df_nonEff[rownames(rowData(snpRSE)),]),  outFile2)

indUnique = which(!duplicated(snpDat$rsID))
snpRSE = snpRSE[indUnique, ]
snpRSE = snpRSE[order(snpRSE)]
snpDat = snpDat[indUnique, ]

###############
# make SNP tracks for pyGenomeTracks
# export.bed(rowRanges(snpRSE), file.path('BED','addiction_gwas_all.bed'))

#lead/independent significant SNP track
snpRSE_lead = snpRSE[rowData(snpRSE)$rsID == rowData(snpRSE)$IndSigSNP,]
rowRanges(snpRSE_lead)$score = -log10(rowRanges(snpRSE_lead)$gwasP)
# export.bed(rowRanges(snpRSE_lead), file.path('BED','addiction_gwas_IndSigSNP.bed'))

#finemapped causal SNPs
snpRSE_fine = snpRSE[rowData(snpRSE)$isCausal,]
rowRanges(snpRSE_fine)$score = rowRanges(snpRSE_fine)$avgPP * 1000 # finemapp PP
# export.bed(rowRanges(snpRSE_fine), file.path('BED','addiction_gwas_finemapped.bed'))

#SNPs in NeuN + peaks
snpRSE_peaks = snpRSE[rowData(snpRSE)$numPeakOverlap > 0,]
rowRanges(snpRSE_peaks)$score = rowRanges(snpRSE_peaks)$numPeakOverlap
# export.bed(rowRanges(snpRSE_peaks), file.path('BED','addiction_gwas_inNeuNPeaks.bed'))

#pre-ML model filtered SNPs
snpRSE_filtered = snpRSE[rowData(snpRSE)$numPeakOverlap > 0 & 
                           rowData(snpRSE)$isCausal,]
length(unique(rowData(snpRSE_filtered)$GenomicLocus)) # 54

####################################################
# inPeak2 SNP in OFC, VLPFC, DLPFC, STC, PUT, NAC
ctx_regions = c('DLPFC_NeuN+','VLPFC_NeuN+','OFC_NeuN+','STC_NeuN+')
str_regions = c('PUT_NeuN+','NAC_NeuN+')

rowData(snpRSE_filtered)$inCtx = apply(assays(snpRSE_filtered)$overlap[,ctx_regions], 1, sum) > 0
rowData(snpRSE_filtered)$inStr = apply(assays(snpRSE_filtered)$overlap[,str_regions], 1, sum) > 0
rowData(snpRSE_filtered)$inPeak2 = rowData(snpRSE_filtered)$inStr | rowData(snpRSE_filtered)$inCtx

####################################################
# define candidate SNPs, w/ any enhancer Z-score > 2
cutOff_effect = .5
cutOff_delta = 0.05

indCtxList = lapply(names(df_effect)[1], function(x) {
  eff = df_effect[rownames(snpRSE_filtered),x]
  non = df_nonEff[rownames(snpRSE_filtered),x]
  # either effect or non-effect allele is predicted enhancer
  isEnh = (eff > cutOff_effect )| (non > cutOff_effect)
  del = df_delta[rownames(snpRSE_filtered),x]
  # SNP is in a cortical NeuN+ peak
  inCtx = rowData(snpRSE_filtered)$inCtx
  # there is a big deltaSNP score
  ind = which(inCtx & isEnh & abs(del) > cutOff_delta)
  return(ind)
  })
names(indCtxList) = names(df_effect)[1]

indStrList = lapply(names(df_effect)[2:3], function(x) {
  eff = df_effect[rownames(snpRSE_filtered),x]
  non = df_nonEff[rownames(snpRSE_filtered),x]
  # either effect or non-effect allele is predicted enhancer
  isEnh = (eff > cutOff_effect )| (non > cutOff_effect)
  del = df_delta[rownames(snpRSE_filtered),x]
  # SNP is in a striatum NeuN+ peak
  inStr = rowData(snpRSE_filtered)$inStr
  # there is a big deltaSNP score
  ind = which(inStr & isEnh & abs(del) > cutOff_delta)
  return(ind)
})
names(indStrList) = names(df_effect)[2:3]

indList = c(indCtxList, indStrList)
lengths(indList)

for(name in names(indList)){
  export.bed(rowRanges(snpRSE_filtered[indList[[name]],]), 
             file.path('BED',paste0('addiction_gwas_candidate_',name,'.bed')))
}

# get all candidate SNPs
snpRSE_ctx_candidate = snpRSE_filtered[unique(unlist(indCtxList)), ]
snpRSE_str_candidate = snpRSE_filtered[unique(unlist(indStrList)), ]
snpRSE_candidate = snpRSE_filtered[unique(unlist(indList)), ]

export.bed(rowRanges(snpRSE_ctx_candidate), file.path('BED','addiction_gwas_ctx_candidate.bed'))
export.bed(rowRanges(snpRSE_str_candidate), file.path('BED','addiction_gwas_str_candidate.bed'))
export.bed(rowRanges(snpRSE_candidate), file.path('BED','addiction_gwas_candidate.bed'))

rsID = unique(unlist(sapply(indList, names)))
write.table(rsID, col.names = F, row.names = F, quote = F,
      file = paste0('tables/candidate_snp_ctx_str_delta',cutOff_delta,'.txt'))


##########################################
# get loci to plot around candidate SNPs
snpBlock1 = GRanges(unique(rowRanges(snpRSE_ctx_candidate)$GenomicLocus))
start(snpBlock1) = floor(start(snpBlock1)/10^3 - 1)*10^3+1 # round to nearest 10kb
end(snpBlock1) = ceiling(end(snpBlock1)/10^3 + 1)*10^3 # round to nearest 10kb
export.bed(snpBlock1, file.path('BED','ctx_candidate_snpBlock.bed'))

# get loci to plot around candidate SNPs
snpBlock2 = GRanges(unique(rowRanges(snpRSE_str_candidate)$GenomicLocus))
start(snpBlock2) = floor(start(snpBlock2)/10^3 - 1)*10^3+1 # round to nearest 10kb
end(snpBlock2) = ceiling(end(snpBlock2)/10^3 + 1)*10^3 # round to nearest 10kb
export.bed(snpBlock2, file.path('BED','str_candidate_snpBlock.bed'))

# get loci to plot around candidate SNPs
snpBlock3 = GRanges(unique(rowRanges(snpRSE_candidate)$GenomicLocus))
start(snpBlock3) = floor(start(snpBlock3)/10^3 - 1)*10^3+1 # round to nearest 10kb
end(snpBlock3) = ceiling(end(snpBlock3)/10^3 + 1)*10^3 # round to nearest 10kb
export.bed(snpBlock3, file.path('BED','snp_candidate_snpBlock.bed'))


###########################################################
# filter GTF file just around snpBlock, for fast plotting
gtf = import('gene_tracks/gencode.v32.basic.annotation.gtf.gz')

# filter out transcript_type and transcript_support_level
# https://www.gencodegenes.org/pages/data_format.html
gtf = gtf[gtf$transcript_support_level %in% c(1, 2)]
gtf = gtf[gtf$transcript_type %in% c('protein_coding', 'retained_intron', 
                                     'lncRNA','polymorphic_pseudogene')]

# get loci to plot around candidate SNPs
snpBlock = GRanges(unique(rowRanges(snpRSE_filtered)$GenomicLocus))
# round to nearest 10kb
start(snpBlock) = floor(start(snpBlock)/5/10^4 - 1)*5*10^4+1 - 10^5 
# round to nearest 10kb
end(snpBlock) = ceiling(end(snpBlock)/5/10^4 + 1)*5*10^4 + 10^5
oo = findOverlaps(query = gtf, subject = snpBlock)
gtf_filtered = gtf[unique(queryHits(oo))]
export(gtf_filtered, 'gene_tracks/gencode.v32.basic.annotation.filtered.gtf')



