library(feather)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(reshape2)
library(SummarizedExperiment)
library(ggsignif)

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

# load the models
celltypes = c('EXC','PV','SST','VIP','D1','D2')
tissue = c('Ctx','Cpu','Nac')
group = c('Ctx-PV','Ctx-SST','Ctx-VIP','Ctx-EXC','Str-PV','Str-SST','Str-D1','Str-D2')

# read in the CNN models
models = read_feather('rdas/SingleTask_IDR_peakPrediction_CNN_performance.feather')
models$sample2 = ss(models$sample, '_fold')
models$celltypes = factor(ss(ss(models$sample,'_',2), 'pos'), celltypes)
models$tissue = factor(ss(models$sample,'_',3), tissue)
models$tissue2 = ifelse(models$tissue== 'Ctx','Ctx','Str')
models$group = factor(paste(models$tissue2, models$celltypes, sep = '-'), group)


# load the model predictions of the SNPs, rotate matrix, and make into numeric
df1 = read.delim('rdas/SingleTask_IDR_addiction_snps_allele_effect_predictions.tsv', 
                 row.names = 1)[-c(1:2),] %>% apply(1, as.numeric)
df2 = read.delim('rdas/SingleTask_IDR_addiction_snps_allele_nonEffect_predictions.tsv', 
                 row.names = 1)[-c(1:2),] %>% apply(1, as.numeric)


# get model names and folds, collapse folds/samples across tissue/celltype
splitInd = split(seq(nrow(models)), models$group)


# average model prediction across folds of same data
df_effect = do.call('cbind',lapply(splitInd, function(i) colMeans(df1[i,]))) %>% data.frame()
df_effect = df_effect[snpDat$rsID, ]
df_nonEff = do.call('cbind',lapply(splitInd, function(i) colMeans(df2[i,]))) %>% data.frame()
df_nonEff = df_nonEff[snpDat$rsID, ]

# compute the delta SNP scores
names(df_effect) = names(df_nonEff) = names(splitInd)
all.equal(rownames(df_effect), rownames(df_nonEff))
df_delta = df_effect - df_nonEff

apply(df_effect, 2, summary)
apply(df_nonEff, 2, summary)
apply(df_delta, 2, summary)

# how many SNPs predicted to be active in
table(apply(df_effect, 1, function(x) sum(x > .5)))
table(apply(df_nonEff, 1, function(x) sum(x > .5)))
table(apply((df_effect > .5) - (df_nonEff > .5), 1, function(x) sum(abs(x) > 0)))



########################################################
# for NeuN+ peak, does fineMapping stratify SNP scores ?
snp_scoreMat_long = melt(df_effect %>% rownames_to_column('rsID'), id.vars = 'rsID',
                         value.name = 'score', variable.name = 'model_type')
snp_scoreMat_long$tissue = ss(as.character(snp_scoreMat_long$model_type),'-',1)
snp_scoreMat_long$celltype = ss(as.character(snp_scoreMat_long$model_type),'-',2)

snp_scores = cbind(as.data.frame(snp_scoreMat_long), 
                   snpDat[match(snp_scoreMat_long$rsID, snpDat$rsID),] )
snp_scores = snp_scores[, !duplicated(names(snp_scores))]
# save(df_effect, df_nonEff, df_delta, snp_scores, file = 'rdas/addiction_snp_cnn_score_20200727.rda')

cortical_NeuN = c('DLPFC_NeuN+','VLPFC_NeuN+','OFC_NeuN+','STC_NeuN+')
cortical_Glia = c('DLPFC_NeuN-','VLPFC_NeuN-','OFC_NeuN-','STC_NeuN-')
cortical_models = c('Ctx-EXC','Ctx-PV','Ctx-SST','Ctx-VIP')
striatal_NeuN = c('PUT_NeuN+','NAC_NeuN+')
striatal_Glia = c('PUT_NeuN-','NAC_NeuN-')
striatal_models = c('Str-PV','Str-SST','Str-D1','Str-D2')

# cortical regions
tmp1 = subset(snp_scores, model_type %in% cortical_models)
tmp1$inNeuN = rowSums(tmp1[,cortical_NeuN]) > 0
tmp1$NeuNneg = rowSums(tmp1[,cortical_Glia]) > 0
# tmp1 = tmp1[tmp1$inNeuN | tmp1$NeuNneg,]

# striatal regions
tmp2 = subset(snp_scores, model_type %in% striatal_models)
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

################################
# plot finemapped & NeuN+ peaks
pdf('plots/snp_scores_peaks_causal_20200727.pdf', width = 5.5, height = 3)
ggplot(tmp, aes(x = inPeak, y = score, fill = inPeak)) + 
  geom_hline(yintercept = 0, color = 'darkgray') + 
  geom_violin(scale='width', draw_quantiles = T) + 
  geom_boxplot(width=0.4) +# coord_flip() + 
  xlab('SNP in a NeuN + or - peak') +
  ylab('ML model score') + 
  facet_wrap( ~ model_type, nrow = 2) + 
  theme_bw(base_size = 11) + 
  scale_fill_manual(values = c('gray', '#F8766D','#00BFC4', 'grey40'), name = 'SNP in peak') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# theme(legend.position = "none")
dev.off()



my_comparisons <- list( c("none", "NeuN-"), c("NeuN+", "NeuN-"), c("none", "NeuN+"), 
                        c("Both", "NeuN+"), c("Both", "NeuN-"))
N = 8*6
pdf('plots/snp_scores_peaks_causal_withPValues_20200730.pdf', width = 5.5, height =3)
ggplot(tmp, aes(x = inPeak, y = score, fill = inPeak)) + 
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
  facet_wrap( ~ model_type, nrow = 2) + 
  scale_fill_manual(name = 'SNP in peak', 
                    values = c('gray', '#F8766D','#00BFC4', 'grey40')) +
  xlab('SNP in a NeuN + or - peak') +
  ylab('ML model score') + 
  ylim(0,1.85) +
  theme_bw(base_size = 10) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


