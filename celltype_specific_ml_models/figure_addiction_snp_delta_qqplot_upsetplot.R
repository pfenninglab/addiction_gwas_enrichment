library(feather)
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggpubr)
library(reshape2)
library(SummarizedExperiment)
library(ggsignif)
library(qqman)
library(qvalue)

options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
groups = c('Ctx-EXC','Cpu-D1','Cpu-D2', 'Nac-D1', 'Nac-D2')

####################################################
# load in annotated SNPs w/ CNN scores for each SNP
out_rds = '../gwas_enrichment/rdas/addiction_snp_cnn_scores_qval_20210413.rds'
df_snp_long = readRDS(out_rds) 
df_snp_long2 = df_snp_long %>% mutate(tmp = paste(group, rsID)) %>%
  filter(!duplicated(tmp)) %>% arrange(snp_delta_pval) %>%
  mutate( group = factor(group, groups),
    observed_log10pval = -log10(snp_delta_pval), 
    expected_log10pval = -log10(ppoints(length(snp_delta_pval))))
alpha = 0.05



########################################
# 1) make qq plot of the delta SNP scores
num.signif  = df_snp_long2 %>% distinct(rsID,group,.keep_all = T) %>% group_by(group) %>% 
  summarize(n=paste0("TierA = ", sum(tier =='TierA'),'\n',
                     "TierB = ", sum(tier =='TierB')))
pdf('figures/addiction_snp_delta_qqplot_20210414.pdf', width = 5.5, height = 2)
ggplot(df_snp_long2,aes(x = expected_log10pval, y = observed_log10pval)) +
  # geom_point(pch = 20, aes(color = factor(snp_delta_qval < alpha))) +
  geom_hex() + scale_fill_viridis_c(trans = "log10") + 
  scale_color_manual(values = c('black', 'red'), name = paste('qval <',alpha)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  xlab(expression(paste('-log'[10],'(Exp ',Delta,' SNP pval)'))) + 
  ylab(expression(paste('-log'[10],'(Obs ',Delta,' SNP pval)'))) + 
  theme_bw(base_size = 8) + facet_grid(~ group) + 
   geom_text(data=num.signif, aes(x=2, y=80, label=n), size = 3, 
            colour="darkred", inherit.aes=FALSE, parse=FALSE) +
  theme(legend.position="right", legend.margin=margin(-10, 0, 0, 0))
dev.off()



pdf('plots/addiction_snp_score_probability_calibration_20210414.pdf', width = 5.5, height = 2)
ggplot(df_snp_long2,aes(x = effect_allele_score, y = effect_allele_prob)) +
  # geom_point(pch = 20, aes(color = factor(snp_delta_qval < alpha))) +
  geom_hex() + scale_fill_viridis_c(trans = "log10") + 
  scale_color_manual(values = c('black', 'red'), name = paste('qval <',alpha)) + 
  xlab('SNP Scores of effect allele') + 
  ylab('Calibrated SNP Pr(active)') + 
  theme_bw(base_size = 8) + facet_grid(~ group) +
  theme(legend.position="right", legend.margin=margin(-10, 0, 0, 0))

ggplot(df_snp_long2,aes(x = snp_delta_prob)) +
  # geom_point(pch = 20, aes(color = factor(snp_delta_qval < alpha))) +
  geom_density() + 
  scale_color_manual(values = c('black', 'red'), name = paste('qval <',alpha)) + 
  xlab('SNP Scores of effect allele') + 
  theme_bw(base_size = 8) + facet_grid(~ group) +
  theme(legend.position="right", legend.margin=margin(-10, 0, 0, 0))
dev.off()




##############################################################
# 2) summary numbers of SNPs across different criteria groups
## number SNPs in NeuN+ peaks
df_snp_long2 %>% filter(inNeuN) %>% distinct(rsID) %>% tally() # 557 snps
df_snp_long2 %>% filter(inNeuN) %>% distinct(GenomicLocus) %>% tally() # 124 loci

## number of qval < 0.05 SNPs and loci
df_snp_long2 %>% filter(tier == 'TierA') %>% distinct(rsID,.keep_all = T) %>% 
  group_by(group) %>% tally()
df_snp_long2 %>% distinct(rsID,.keep_all = T) %>%  group_by(tier) %>% tally() # 55, 505, 502
df_snp_long2 %>% filter(tier == 'TierA') %>% distinct(GenomicLocus) %>% tally() # 37 loci
df_snp_long2 %>% filter(tier == 'TierB') %>% distinct(GenomicLocus) %>% tally() # 118 loci
df_snp_long2 %>% filter(tier == 'TierC') %>% distinct(GenomicLocus) %>% tally() # 121 loci

## 37 loci w/ NeuN+ signif SNPs
df_snp_long2 %>% filter(snp_delta_qval < alpha) %>%
  mutate(inNeuN = rowSums(across(ends_with("NeuN+"))) > 0) %>% 
  filter(inNeuN) %>% distinct(GenomicLocus) %>% tally()

## 55 signif SNPs in NeuN+ peaks
df_snp_long2 %>% filter(snp_delta_qval < alpha) %>%
  mutate(inNeuN = rowSums(across(ends_with("NeuN+"))) > 0) %>% 
  filter(inNeuN) %>% distinct(rsID) %>% tally()

## number of qval < 0.05 SNPs and loci
df_snp_long2 %>% filter(snp_delta_qval < alpha) %>%
  mutate(inNeuN = rowSums(across(ends_with("NeuN+"))) > 0) %>% 
  filter(inNeuN) %>% group_by(group) %>% tally()

#####################################################################
# make upset plot of predicted delta SNPs by model and NeuN+ overlap
library(ggupset)
# SNPs overlapping NeuN+ peaks
tmp1 = df_snp_long2 %>% filter(snp_delta_qval < alpha) %>%
  mutate(inNeuN = rowSums(across(ends_with("NeuN+"))) > 0, 
         group = 'inNeuN+') %>% filter(inNeuN) %>% 
  distinct(rsID, .keep_all=TRUE) %>% select(c(rsID, group))

# SNPs w/ significant delta SNP score
tmp2 = df_snp_long2 %>% filter(snp_delta_qval < alpha) %>%
  mutate(group = as.character(group)) %>%
  select(c(rsID, group)) 

pdf('figures/addiction_snp_delta_upset_20210413.pdf', width = 11, height = 2.5)
rbind(tmp1, tmp2) %>% mutate(group = factor(group, c(groups, 'inNeuN+'))) %>% 
  arrange(group) %>% group_by(rsID) %>% summarise(group = list(group)) %>%
  ggplot(aes(x=group)) + geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,size = 2) +
  scale_x_upset(order_by = 'degree') +
  scale_y_continuous(breaks = NULL, name = "#SNPs qval < 0.05", 
                     lim = c(0, 200)) + 
  theme_bw(base_size = 8) +
  theme_combmatrix(combmatrix.label.text = element_text(size=8),
                   combmatrix.label.extra_spacing = 5)
dev.off()

# SNPs w/ significant delta SNP score
tmp2 = df_snp_long2 %>% filter(snp_delta_qval < alpha) %>%
  mutate(group = as.character(group)) %>%
  select(c(rsID, group)) 

pdf('figures/addiction_snp_delta_upset_20210413.pdf', width = 11, height = 2.5)
rbind(tmp1, tmp2) %>% mutate(group = factor(group, c(groups, 'inNeuN+'))) %>% 
  arrange(group) %>% group_by(rsID) %>% summarise(group = list(group)) %>%
  ggplot(aes(x=group)) + geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,size = 2) +
  scale_x_upset(order_by = 'degree') +
  scale_y_continuous(breaks = NULL, name = "#SNPs qval < 0.05", 
                     lim = c(0, 200)) + 
  theme_bw(base_size = 8) +
  theme_combmatrix(combmatrix.label.text = element_text(size=8),
                   combmatrix.label.extra_spacing = 5)
dev.off()

