library(feather)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(reshape2)
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#########################
# organize the celltype
# celltype = c('EXC','PV','VIP','SST','D1','D2')
# colors = c('#e41a1c','#ff7f00','#984ea3', '#ffff33','#377eb8','#4daf4a')
celltype = c('EXC','D1','D2')
colors = c('#e41a1c','#377eb8','#4daf4a')
names(colors) = celltype
tissue = c("Ctx", "Cpu", "Nac")
lab = c("Mo2015", "Pfenning")

####################
# read in the data
df = read.delim('tables/SingleTask_IDR_peakPrediction_CNN_testSet_performance.tsv')
df$lab = factor(ss(basename(df$model), '_',1), lab)
df$celltype = factor(gsub('pos','',ss(basename(df$model), '_',2)), celltype)
df$tissue = factor(ss(basename(df$model), '_',3), tissue)

pdf('plots/mouse_ml_model-CNN_testSet_performance.pdf', width = 11, height = 5)
acc = ggplot(data = df, aes(x = celltype, y = accuracy)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_y_continuous(limits = c(.5, NA)) + # baseline accuracy is 50%
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('Accuracy')

f1s = ggplot(data = df, aes(x = celltype, y = f1_score)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_y_continuous(limits = c(0, NA)) + # baseline f-score is 0
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('F1 score')


roc = ggplot(data = df, aes(x = celltype, y = auROC)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_y_continuous(limits = c(.5, NA)) + # baseline auROC is 50%
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('auROC')


prc = ggplot(data = df, aes(x = celltype, y = auPRC)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_y_continuous(limits = c(0, NA)) + # baseline f-score is 0
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('auPRC')

fhs = ggplot(data = df, aes(x = celltype, y = fhalf_score)) +
  geom_boxplot(aes(fill = celltype)) +
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) +
  scale_fill_manual(values = colors, labels = celltype) +
  scale_y_continuous(limits = c(0, NA)) + # baseline f-score is 0
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') +
  theme_bw(base_size = 10) + ggtitle('F-half score')

acc; prc; roc; f1s; fhs
dev.off()


#################################
# make performance figures
df = df[complete.cases(df), ]
df_long = melt(df, id.vars =c('sample', 'index', 'model_class','model',
                              'model_type','model_species','lab','celltype','tissue'))
df_long = df_long[df_long$variable %in% c('accuracy','auPRC','auROC','f1_score'),]

pdf('figures/figure_mouse_ml_model-CNN_testSet_performance.pdf', 
    width = 5.5, height = 5)
ggplot(data = df_long, aes(x = celltype, y = value)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(width = .25, pch = 21, aes(fill = celltype)) + 
  ylim(c(0.5,NA)) + 
  scale_fill_manual(values = colors, labels = celltype) + 
  facet_grid(variable~lab + tissue, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = 11)
dev.off()

