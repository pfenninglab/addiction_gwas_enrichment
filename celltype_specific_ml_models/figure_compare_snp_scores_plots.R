options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(RColorBrewer)
library(ggplot2)
library(feather)
library(reticulate)
library(umap)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)

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

############################################
# load in CNN snp scores for D1 D2 CNN predictions
cnn_scores = feather::read_feather('rdas/addiction_snps_D1-D2-PV-EXC_VIP_CNN_scores.feather')
cnn_scores$allele = toupper(ss(cnn_scores$SNP,'_',2))
cnn_scores$SNP = ss(cnn_scores$SNP,'_',1)

#############################
# read in gkmSVM scores
gkmSVM_files = list.files(path = 'gkmSVM',pattern = 'addiction_snps_gkmSVM_', full.names = T)
gkm_scores = do.call('bind_rows',lapply(gkmSVM_files, function(x){
  tmp = tibble:: as_tibble(read.delim(x, col.names = c('SNP','score'), header = F))
  tmp$sample = 'Addiction_snps'
  tmp$model_class = 'SVM'
  tmp$model_type = ss(ss(x, '_',5),'\\.')
  tmp$model_species = 'Mouse'
  tmp$model = x
  tmp$test_species = 'Human'
  tmp$allele = ss(x, '_',4)
  return(tmp)
}))

######################
# combine snp scores 
snp_scores = bind_rows(cnn_scores, gkm_scores)
snp_scores$model_type[snp_scores$model_type=='PVPOSvsPVNEG'] = 'PVvsPVNEG'
snp_scores2 = snp_scores %>% 
  group_by(sample, model_class, model_type, model_species, SNP ) %>%
  summarize(score = mean(score))

#############################################################
# reshape to have snp by model prediction & compute Z scores
snp_scoreMat = acast(snp_scores2,SNP ~ model_class + model_type , value.var = 'score')
snp_ZscoreMat = apply(snp_scoreMat, 2, function(x) (x-mean(x))/sd(x)) # znormalize by model
snp_ZscoreMat_long = melt(snp_ZscoreMat, varnames = c('SNP', 'model'), value.name = 'zscore')
snp_ZscoreMat_long$model_class = ss(as.character(snp_ZscoreMat_long$model), '_', 1)
snp_ZscoreMat_long$model_type = ss(as.character(snp_ZscoreMat_long$model), '_', 2)

########################################################
# for NeuN+ peak, does fineMapping stratify SNP scores ?
snp_scores3 = cbind(as.data.frame(snp_ZscoreMat_long), 
                    snpDat[match(snp_ZscoreMat_long$SNP, snpDat$rsID),] )
snp_scores3 = snp_scores3[, !duplicated(names(snp_scores3))]
# save(snp_ZscoreMat, snp_scores3, file = 'rdas/addiction_snp_cnn_svm_zscore.rda')

cortical_NeuN = c('DLPFC_NeuN+','VLPFC_NeuN+','OFC_NeuN+','STC_NeuN+')
cortical_Glia = c('DLPFC_NeuN-','VLPFC_NeuN-','OFC_NeuN-','STC_NeuN-')
cortical_models = c('PVvsPVNEG', 'PVvsEXC')
striatal_NeuN = c('PUT_NeuN+','NAC_NeuN+')
striatal_Glia = c('PUT_NeuN-','NAC_NeuN-')
striatal_models = c('D1D2vsBulk','D1vsBulk','D2vsBulk')

# cortical regions
tmp1 = subset(snp_scores3, model_type %in% cortical_models)
tmp1$inNeuN = rowSums(tmp1[,cortical_NeuN]) > 0
tmp1$NeuNneg = rowSums(tmp1[,cortical_Glia]) > 0
tmp1 = tmp1[tmp1$inNeuN | tmp1$NeuNneg,]

# striatal regions
tmp2 = subset(snp_scores3, model_type %in% striatal_models)
tmp2$inNeuN = rowSums(tmp2[,striatal_NeuN]) > 0
tmp2$NeuNneg = rowSums(tmp2[,striatal_Glia]) > 0
tmp2 = tmp2[tmp2$inNeuN | tmp2$NeuNneg,]

tmp = rbind(tmp1, tmp2)
tmp$isCausal = factor(ifelse(tmp$isCausal, 'isCausal','notCausal'),
                      levels = c('notCausal','isCausal'))
tmp$inPeak = factor(ifelse(tmp$inNeuN & tmp$NeuNneg, 'inBoth',
             ifelse(tmp$inNeuN & !tmp$NeuNneg, 'inNeuN', 'inGlia')), 
             levels = c('inGlia','inBoth','inNeuN'))
tmp = subset(tmp, tmp$inPeak %in% c('NeuN-','NeuN+'))

################################
# plot finemapped & NeuN+ peaks
pdf('plots/snp_zscores_peaks_causal_20200707.pdf', width = 5.5, height = 3)
ggplot(tmp, aes(x = inPeak, y = zscore, fill = inPeak)) + 
  geom_hline(yintercept = 0, color = 'darkgray') + 
  geom_violin() + geom_boxplot(width=0.4) +# coord_flip() + 
  xlab('SNP in NeuN+/- only peak') +
  ylab('ML model zscore of SNP') + 
  facet_grid(model_class ~ model_type, scales = 'free') + 
  theme_bw(base_size = 11) + 
  scale_fill_discrete(name = '') + 
  theme(legend.position = "none")
dev.off()


#########################################################################
### linear model of shift in model predictions by being in peak & causal
tmpList = split(tmp, factor(tmp$model))
lmList = lapply(tmpList, function(x){
  mod = lm(zscore ~ inPeak, data = x)
  summary(mod)
})

# transform linear model outputs to data frame
lmDat = do.call('rbind', lapply(lmList,function(x) data.frame( x$coefficients)))
lmDat = lmDat %>% rownames_to_column() %>%
  # split rowname into ML model & test
  mutate(model =  ss(rowname, '\\.', 1), 
         effect = ss(rowname, '\\.', 2),
         model_class = ss(as.character(model), '_', 1),
         model_type = ss(as.character(model), '_', 2)) %>%
  filter(effect != "(Intercept)") %>%
  group_by(model, effect) %>% 
  mutate(FDR = p.adjust(`Pr...t..`, 'BH'),
         FDR_signif = FDR < 0.05)

lmDat$effect = gsub('inPeak','',lmDat$effect )
lmDat$effect = gsub('isCausalisCausal','isCausal',lmDat$effect )
lmDat$effect = factor(lmDat$effect, levels = rev(c('inNeuN','isCausal','inNeuN:isCausal')))

lmDat %>% group_by(model_class, model_type) %>%
  filter(FDR < 0.01) %>% summarise(FDR)

####################################
### plot results of linear models
pdf('plots/snp_lmod_peaks_causal_20200707.pdf', width = 5.5, height = 2.5)
ggplot(lmDat, aes(x = effect, y = Estimate, alpha = FDR_signif)) + 
  geom_bar(stat = 'identity') + coord_flip() + 
  geom_errorbar(aes(ymin= (Estimate-Std..Error), ymax=(Estimate+Std..Error)), 
                width = .3, position=position_dodge(width=.9)) +
  scale_color_manual(values=c(NA,'black')) +
  scale_alpha_manual(values=c(.3,1), name = 'FDR < 0.05')+
  guides(color = 'none') +
  geom_hline(yintercept = 0, colour = 'black') + 
  facet_grid( model_class ~ model_type) + 
  ylab('Effect size (ML model zscore)') +
  xlab('Linear Model Effect') +
  theme_bw(base_size = 10) + 
  theme(axis.text.y = element_text(face="bold", size=6), 
        axis.text.x = element_text(size=6), 
        legend.text=element_text(size=6),
        legend.title=element_text(size=6))
dev.off()





