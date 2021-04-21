library(feather)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(reshape2)
library(wesanderson)
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#########################
# organize the celltype
celltype = c('EXC','PV','VIP','SST','D1','D2')
colors = c('#e41a1c','#ff7f00','#984ea3', '#ffff33','#377eb8','#4daf4a')
names(colors) = celltype
tissue = c("Ctx", "Cpu", "Nac")
lab = c("Mo2015", "Pfenning")

####################
# read in the data
tmp1 = read_tsv('tables/SingleTask_IDR_peakPrediction_CNN_v2_performance.tsv' , 
           col_type = cols()) %>% mutate(version = 'v2') %>% select(-X1)
  
tmp2 = read_tsv('tables/SingleTask_IDR_peakPrediction_CNN_v3_performance.tsv', 
               col_type = cols()) %>% mutate(version = 'v3') %>% select(-X1)
df = rbind(tmp1,tmp2)
df <- df[,colSums(is.na(df))<nrow(df)]
df$lab = factor(ss(basename(df$model), '_',1), lab)
df$celltype = factor(gsub('pos','',ss(basename(df$model), '_',2)), celltype)
df$tissue = factor(ss(basename(df$model), '_',3), tissue)
df$dropout = factor(df$dropout)

pdf('plots/mouse_ml_model-CNN_v3_performance.pdf', width = 11, height = 5)
acc = ggplot(data = df, aes(x = celltype, y = accuracy, fill = dropout)) + 
  geom_boxplot(aes(color = celltype),lwd=.5) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(color = celltype)) + 
  scale_color_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_fill_nejm()+ 
  scale_y_continuous(limits = c(.5, NA)) + # baseline accuracy is 50%
  facet_grid(version~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('Accuracy')

f1s = ggplot(data = df, aes(x = celltype, y = f1_score, fill = dropout)) + 
  geom_boxplot(aes(color = celltype),lwd=.5) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(color = celltype)) + 
  scale_color_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_fill_nejm()+ 
  scale_y_continuous(limits = c(0.5, NA)) + # baseline f-score is 0
  facet_grid(version~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('F1 score')


roc = ggplot(data = df, aes(x = celltype, y = auROC, fill = dropout)) + 
  geom_boxplot(aes(color = celltype),lwd=.5) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(color = celltype)) + 
  scale_color_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_fill_nejm()+ 
  scale_y_continuous(limits = c(.5, NA)) + # baseline auROC is 50%
  facet_grid(version~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('auROC')


prc = ggplot(data = df, aes(x = celltype, y = auPRC, fill = dropout)) + 
  geom_boxplot(aes(color = celltype),lwd=.5) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(color = celltype)) + 
  scale_color_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_fill_nejm()+ 
  scale_y_continuous(limits = c(0.5, NA)) + # baseline f-score is 0
  facet_grid(version~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('auPRC')

fhs = ggplot(data = df, aes(x = celltype, y = fhalf_score, fill = dropout)) +
  geom_boxplot(aes(color = celltype),lwd=.5) +
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(color = celltype)) +
  scale_color_manual(values = colors, labels = celltype, guide = 'none') + 
  scale_fill_nejm()+ 
  scale_y_continuous(limits = c(0.5, NA)) + # baseline f-score is 0
  facet_grid(version~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('F-half score')

acc; prc; roc; f1s; fhs
dev.off()

#################################
# make performance figures
celltypes = c('EXC','D1','D2')
colors = c('#e41a1c','#377eb8','#4daf4a')
df2 = df %>% filter(version == 'v3') %>% filter(dropout == 0.25) %>%
  mutate(celltype = factor(gsub('pos','',ss(basename(model), '_',2)), celltypes)) %>% 
  filter(!is.na(celltype))

df_long = melt(df2, id.vars =c('sample', 'index','model',
                              'model_type','model_species','lab','celltype','tissue'))
df_long = df_long[df_long$variable %in% c('accuracy','auPRC','auROC','f1_score'),]
df_long$value = as.numeric(df_long$value)

out_fn = 'tables/SingleTask_IDR_bestModels_CNN_performance.tsv'
write.table(df2, file = out_fn, quote = F, sep = '\t', row.names = F)

pdf('figures/figure_mouse_ml_model-CNN_v3_performance.pdf', width = 5.5, height = 2)
ggplot(data = df2, aes(x = celltype, y = auPRC, fill = celltype)) + 
  geom_boxplot(aes(fill = celltype),lwd=.5) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltype) + 
  scale_y_continuous(limits = c(0.5, NA)) + # baseline PRC is fraction positives
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10)
dev.off()

