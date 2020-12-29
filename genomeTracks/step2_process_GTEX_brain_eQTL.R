options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyverse)

######################
# load in annotated SNPs
load('../gwas_enrichment/rdas/addiction_snp_withCausalDB_Fullard_overlap_rse.rda')
indKeep = with(snpDat, which(!duplicated(rsID) & isCausal))
indKeep = with(snpDat, which(!duplicated(rsID) & isCausal & numPeakOverlap>0))
snpDat = snpDat[indKeep, ]
snpDat$gtex_id = with(snpDat, paste(paste0('chr',chr), pos, sep = '_'))

############################################
# read in variants w/ signif eQTL relation
sigVarFiles = list.files(path = 'GTEX_Analysis_v8_eQTL', 
                   pattern='signif_variant_gene_pairs.txt.gz',full.names = T)
names(sigVarFiles) = gsub('Brain_', '',ss(basename(sigVarFiles), '\\.'))
sigVarList = lapply(sigVarFiles, read_tsv)
sigVarMat = bind_rows( sigVarList, .id = "Region")

# add back chromosome position columns
sigVarMat$chr = ss(sigVarMat$variant_id, '_', 1)
sigVarMat$pos = as.numeric(ss(sigVarMat$variant_id, '_', 2))
sigVarMat$gtex_id = with(sigVarMat, paste(chr, pos, sep = '_'))

# match eQTLs to addiction SNPs in NeuN+ peaks and predicted causal
indEQTL = which(sigVarMat$gtex_id %in% snpDat$gtex_id)
sigVarMat = sigVarMat[indEQTL, ]

# add back allele columns
sigVarMat$ref = ss(sigVarMat$variant_id, '_', 3)
sigVarMat$alt = ss(sigVarMat$variant_id, '_', 4)
length(unique(sigVarMat$variant_id))
save(sigVarMat, file = 'rdas/Brain_GTEX_eQTL_filtered_NeuN+_Causal_SNP.rda')

############################################
# make into pyGenomeTrack links formats
max(-log10(sigVarMat$pval_beta)) # 54
linkCTX = with(sigVarMat[sigVarMat$Region %in% c('Anterior_cingulate_cortex_BA24',
                                               'Cortex','Frontal_Cortex_BA9'),], 
             data.frame(chr1 = chr, start1 = pos, end1 = pos + 1, 
                        chr2 = chr, start2 = pos - tss_distance, 
                        end2 = pos -tss_distance+ 1, score = -log10(pval_beta)))
write.table(linkCTX, file = 'BED/GTEX_brain_NeuN+_Causal_eQTL_CTX.txt',
            quote = FALSE, sep = '\t',col.names = F,row.names = F)

linkSTR = with(sigVarMat[sigVarMat$Region %in% c('Caudate_basal_ganglia',
                                                 'Nucleus_accumbens_basal_ganglia',
                                                 'Putamen_basal_ganglia'),], 
               data.frame(chr1 = chr, start1 = pos, end1 = pos + 1, 
                          chr2 = chr, start2 = pos - tss_distance, 
                          end2 = pos-tss_distance+ 1, score = -log10(pval_beta)))
write.table(linkSTR, file = 'BED/GTEX_brain_NeuN+_Causal_eQTL_STR.txt',
            quote = FALSE, sep = '\t',col.names = F,row.names = F)


eQTL_gene_loc = with(sigVarMat, 
               data.frame(chr = chr, start = pos - tss_distance, 
                          end = pos-tss_distance+ 1))
eQTL_gene_loc = eQTL_gene_loc[!duplicated(eQTL_gene_loc), ]
write.table(eQTL_gene_loc, file = 'BED/GTEX_brain_NeuN+_Causal_eQTL_gene_loc.bed',
            quote = FALSE, sep = '\t',col.names = F,row.names = F)

