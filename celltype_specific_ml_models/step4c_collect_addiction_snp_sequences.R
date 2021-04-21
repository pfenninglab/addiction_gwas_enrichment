library(feather)
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggpubr)
library(reshape2)
library(SummarizedExperiment)
library(ggsignif)
library(qvalue)
library(swfdr)

options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

######################
# load in annotated SNPs
load('../gwas_enrichment/rdas/addiction_snp_withCausalDB_Fullard_NeuN+-_overlap_rse.rda')
rowData(snpRSE) = snpDat2 = snpDat
indUnique = which(!duplicated(snpDat$rsID))
snpRSE = snpRSE[indUnique, ]
snpDat = snpDat[indUnique, ]
rownames(snpDat) = snpDat$rsID
table(snpDat$isCausal)
table(snpDat$numPeakOverlap>0)
table(isCausal = snpDat$isCausal, inPeak = snpDat$numPeakOverlap>0)

##################################
# get the 500bp window GC content
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
genome <- BSgenome.Hsapiens.UCSC.hg38.masked
offset = 500/2
seqEffect <-getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos+offset)
snpDat$GC_content = letterFrequency(seqEffect, letters = "GC", as.prob = TRUE)


##############################################
# load the models and prepare metadata
tissues = c('Ctx','Cpu','Nac')
celltypes = c('EXC','D1','D2')
groups = c('Ctx-EXC','Cpu-D1','Cpu-D2', 'Nac-D1', 'Nac-D2')

# read in the CNN models
models = read_tsv('tables/SingleTask_IDR_bestModels_CNN_performance.tsv') %>% 
  mutate(group = factor(paste(tissue, celltype, sep = '-'), groups)) %>% select(-index)

# read the GC content of each model
models = models %>% mutate(
  avgGC = map_dbl(prefix,function(prefix){
    fasta_fn = file.path('FASTA_CV',paste0(prefix,'_validPos.fa'))
    seq = readDNAStringSet(fasta_fn, format="fasta")
    ret = mean(letterFrequency(seq, letters = "GC", as.prob = TRUE))
    return(ret)
  } ))
models2 = models %>% group_by(group) %>% summarise(avgGC = mean(avgGC))

# read the calibration of CNN scores to percent positive
calibration_fn = 'rdas/SingleTask_IDR_CNN_validation_pos_calibration_ecdf.rds'
model_calibration = readRDS(calibration_fn)

#########################################
# load the model predictions of the SNPs
df_effect = read_tsv('tables/SingleTask_IDR_bestModels_addiction_snps_allele_effect_predictions.tsv.gz') %>%
  right_join(models, by = 'model') %>% 
  pivot_wider(id_cols = c('label', 'group'), values_from = 'y_pred_score', 
              names_from ='group', values_fn = mean) %>% 
  column_to_rownames(var = "label") 

df_nonEff = read_tsv('tables/SingleTask_IDR_bestModels_addiction_snps_allele_nonEffect_predictions.tsv.gz') %>%
  right_join(models, by = 'model')  %>% 
  pivot_wider(id_cols = c('label', 'group'), values_from = 'y_pred_score', 
              names_from ='group', values_fn = mean) %>% 
  column_to_rownames(var = "label")

#########################################
# combine SNP predictions across models 
df_effect_long = df_effect %>% rownames_to_column('rsID') %>% 
  pivot_longer(!rsID, names_to = 'group', values_to = 'effect_allele_score') %>%
  group_by(group) %>%
  mutate(effect_allele_prob = model_calibration[[unique(group)]](effect_allele_score)) %>%
  ungroup() %>% mutate(group = factor(group, groups))

df_nonEff_long = df_nonEff %>% rownames_to_column('rsID') %>% 
  pivot_longer(!rsID, names_to = 'group', values_to = 'nonEff_allele_score')%>%
  group_by(group) %>%
  mutate(nonEff_allele_prob = model_calibration[[unique(group)]](nonEff_allele_score)) %>%
  ungroup() %>% mutate(group = factor(group, groups))

toJoin = c('rsID', 'group')
snpDatCols = c('rsID', 'GC_content',"isCausal","numPeakOverlap",'MAF')
df_snp_long = right_join(df_effect_long, df_nonEff_long, by = toJoin) %>%
  mutate(
    snp_delta_score = effect_allele_score - nonEff_allele_score,
    snp_delta_prob = effect_allele_prob - nonEff_allele_prob) %>%
  right_join(snpDat[,snpDatCols], by = 'rsID') %>% 
  right_join(models2, by = "group") %>%
  mutate(diffGC = abs(GC_content - avgGC)) %>% group_by(group) %>%
  # compute the p-value against normal distribution
  mutate(snp_delta_pval = 2* pnorm(q=abs(snp_delta_prob), mean = mean(snp_delta_prob), 
                                   sd = sd(snp_delta_prob), lower.tail=FALSE)) %>%
  ungroup()

# compute the delta SNP scores
df_delta = df_snp_long %>% 
  pivot_wider(id_cols= 'rsID',values_from = 'snp_delta_prob', names_from = 'group')
df_delta = df_delta[rownames(df_nonEff),]

############################################
# compute covariate conditioned FDR w/ swfdr
covariates = c("diffGC", "isCausal","numPeakOverlap",'MAF')
set.seed(1)
lm_qvalue <- lm_qvalue(df_snp_long$snp_delta_pval, X=df_snp_long[,covariates])
lm_qvalue
df_snp_long$snp_delta_qval = lm_qvalue$qvalues

################################################
## merge w/ the SNP across multiple GWAS traits
df_snp_long = df_snp_long %>% left_join(snpDat2) %>% arrange(snp_delta_qval)

## add in tiers
alpha = 0.05
df_snp_long = df_snp_long %>% 
  mutate(inNeuN = rowSums(across(ends_with("NeuN+"))) > 0,
         tier = case_when(inNeuN & snp_delta_qval < alpha ~ 'TierA',
                          snp_delta_qval < alpha ~ 'TierB',
                          inNeuN ~ 'TierC', TRUE ~ 'Other'), 
         tier = factor(tier, c('TierA','TierB','TierC','Other'))) %>%
  relocate(tier, .after = snp_delta_qval)
table(df_snp_long$tier, df_snp_long$group)

out_rds = '../gwas_enrichment/rdas/addiction_snp_cnn_scores_qval_20210413.rds'
saveRDS(df_snp_long, out_rds)
save(df_effect, df_nonEff, df_delta, df_snp_long,
     file = 'rdas/addiction_snp_cnn_score_20210408.rda')
out_xlsx = '../gwas_enrichment/tables/addiction_snp_cnn_scores_qval_20210413.xlsx'
writexl::write_xlsx(df_snp_long, out_xlsx)


###############################################
## export the top SNPs for model interpretation
seqEffect = readDNAStringSet('FASTA/addiction_snps_allele_effect_501.fa', format="fasta")
seqNonEffect = readDNAStringSet('FASTA/addiction_snps_allele_nonEffect_501.fa', format="fasta")

top_snps = df_snp_long %>% filter(tier == 'TierA') %>% filter(!duplicated(rsID))
  
seqEffect = seqEffect[names(seqEffect) %in% top_snps$rsID]
seqNonEffect = seqNonEffect[names(seqNonEffect) %in% top_snps$rsID]

writeXStringSet(seqEffect, paste0('FASTA/top_addiction_snps_allele_effect_501.fa'),
                format = 'fasta',width = 501)
writeXStringSet(seqNonEffect, paste0('FASTA/top_addiction_snps_allele_nonEffect_501.fa'),
                format = 'fasta',width = 501)



########################################################
# for NeuN+ peak, does fineMapping stratify SNP scores ?
# snp_scoreMat_long = df_snp_long %>% 
#   pivot_wider(id_cols = 'rsID', names_from = 'group', values_from = 'y_pred_prob')
# 
# snp_scoreMat_long = melt(df_effect %>% rownames_to_column('rsID'), id.vars = 'rsID',
#                          value.name = 'score', variable.name = 'model_type')
# snp_scoreMat_long$tissue = ss(as.character(snp_scoreMat_long$group),'-',1)
# snp_scoreMat_long$celltype = ss(as.character(snp_scoreMat_long$group),'-',2)
# 
# snp_scores = cbind(as.data.frame(snp_scoreMat_long), 
#                    snpDat[match(snp_scoreMat_long$rsID, snpDat$rsID),] )
# snp_scores = snp_scores[, !duplicated(names(snp_scores))]

cortical_NeuN = c('DLPFC_NeuN+','VLPFC_NeuN+','OFC_NeuN+','STC_NeuN+')
cortical_Glia = c('DLPFC_NeuN-','VLPFC_NeuN-','OFC_NeuN-','STC_NeuN-')
cortical_models = c('Ctx-EXC')
striatal_NeuN = c('PUT_NeuN+','NAC_NeuN+')
striatal_Glia = c('PUT_NeuN-','NAC_NeuN-')
striatal_models = c('Cpu-D1','Cpu-D2', 'Nac-D1', 'Nac-D2')

# cortical regions
tmp1 = subset(df_snp_long, group %in% cortical_models)
tmp1$inNeuN = rowSums(tmp1[,cortical_NeuN]) > 0
tmp1$NeuNneg = rowSums(tmp1[,cortical_Glia]) > 0
# tmp1 = tmp1[tmp1$inNeuN | tmp1$NeuNneg,]

# striatal regions
tmp2 = subset(df_snp_long, group  %in% striatal_models)
tmp2$inNeuN = rowSums(tmp2[,striatal_NeuN]) > 0
tmp2$NeuNneg = rowSums(tmp2[,striatal_Glia]) > 0
# tmp2 = tmp2[tmp2$inNeuN | tmp2$NeuNneg,]

tmp = rbind(tmp1, tmp2)
tmp$isCausal = factor(ifelse(tmp$isCausal, 'isCausal','notCausal'),
                      levels = c('notCausal','isCausal'))
tmp$inPeak = factor(ifelse(tmp$inNeuN & tmp$NeuNneg, 'Both',
                           ifelse(tmp$inNeuN & !tmp$NeuNneg, 'NeuN+', 
                                  ifelse(!tmp$inNeuN & tmp$NeuNneg, 'NeuN-','none'))), 
                    levels = c('none','NeuN-','NeuN+','Both'))
tmp$group = factor(tmp$group, group)

################################
my_comparisons <- list( c("none", "NeuN-"), c("NeuN+", "NeuN-"), c("none", "NeuN+"), 
                        c("Both", "NeuN+"), c("Both", "NeuN-"))
N = 6*3 # 6 possible pairwise comparisons, 3 models

pdf('figures/snp_scores_peaks_causal_withPValues_20210414.pdf', width = 5.5, height =2.5)
ggplot(tmp, aes(x = inPeak, y = pmax(effect_allele_prob, nonEff_allele_prob), fill = inPeak)) + 
  geom_hline(yintercept = 0, color = 'darkgray') + 
  geom_violin(scale='width', draw_quantiles = T) + 
  geom_boxplot(width=0.4) +# coord_flip() + 
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, 
                     p.adjust.method = "bonferroni", group.by = "dose", 
                     label = "p.signif", size = 2.5, 
                     tip.length = .05,
                     bracket.size = 0.3,
                     step.increase = .18,
                     # manual Bonferroni cutoffs
                     symnum.args = list(cutpoints = c(0, 0.001/N, 0.01/N, 0.05/N, 1),
                                        symbols = c("***", "**", "*", "ns")))+ 
  # other plot labels/aesthetics
  facet_wrap( ~ factor(group), nrow = 1) + 
  scale_fill_manual(name = 'SNP in peak', 
                    values = c('gray', '#F8766D','#00BFC4', 'grey40')) +
  xlab('SNP in a NeuN + or - peak') +
  ylab('Calibrated Pr(active CRE)') + 
  ylim(0,1.85) +
  theme_bw(base_size = 10) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

