library(arrow)
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
celltypes = c('EXC','PV','VIP','SST','D1','D2')
colors = c('#e41a1c','#ff7f00','#984ea3', '#ffff33','#377eb8','#4daf4a')
names(colors) = celltypes
labs = c("Mo2015", "Pfenning")
tissues = c("Ctx", "Cpu", "Nac")
groups = c('Ctx-EXC','Cpu-D1','Cpu-D2', 'Nac-D1', 'Nac-D2')

####################
# read in the test performances 
models = read.delim('tables/SingleTask_IDR_bestModels_CNN_performance.tsv')  %>% 
  mutate(group = factor(paste(tissue, celltype, sep = '-'), groups),
         lab = factor(lab, labs), celltype = factor(celltype, celltypes),
         tissue = factor(tissue, tissues))

perf_fn = list.files(path = 'test_performance', pattern = '.feather', recursive = T, full.names = T)
names(perf_fn) = basename(perf_fn)
df = perf_fn %>% lapply(read_feather) %>% bind_rows() %>% 
  right_join(y = models, by ='model_name', suffix = c("", ".validation")) 


fpr_fn = list.files(path = 'predictions', pattern = '.test_negatives.predictions.txt', recursive = T, full.names = T)
names(fpr_fn) = basename(fpr_fn)
fpr = fpr_fn %>% lapply(read_tsv, col_types = cols()) %>% 
  data.table::rbindlist(idcol = 'file') %>% 
  mutate(prefix = ss(file,'_OCP')) %>% 
  group_by(prefix) %>% summarise(fpr = sum(y_pred_score > 0.5) / n())

df = df %>% right_join(fpr)

pdf('plots/mouse_ml_model-CNN_testSet_performance.pdf', width = 11, height = 5)
acc = ggplot(data = df, aes(x = celltype, y = accuracy)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltypes, guide = 'none') + 
  scale_y_continuous(limits = c(.5, NA)) + # baseline accuracy is 50%
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('Accuracy')

f1s = ggplot(data = df, aes(x = celltype, y = f1_score)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltypes, guide = 'none') + 
  scale_y_continuous(limits = c(0, NA)) + # baseline f-score is 0
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('F1 score')


roc = ggplot(data = df, aes(x = celltype, y = auROC)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltypes, guide = 'none') + 
  scale_y_continuous(limits = c(.5, NA)) + # baseline auROC is 50%
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('auROC')


prc = ggplot(data = df, aes(x = celltype, y = auPRC)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = colors, labels = celltypes, guide = 'none') + 
  scale_y_continuous(limits = c(0, NA)) + # baseline f-score is 0
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') + 
  theme_bw(base_size = 10) + ggtitle('auPRC')

fhs = ggplot(data = df, aes(x = celltype, y = fhalf_score)) +
  geom_boxplot(aes(fill = celltype)) +
  geom_jitter(position = position_jitterdodge(), pch = 21, aes(fill = celltype)) +
  scale_fill_manual(values = colors, labels = celltypes) +
  scale_y_continuous(limits = c(0, NA)) + # baseline f-score is 0
  facet_grid(~lab + tissue, scales = 'free_x',space = 'free_x') +
  theme_bw(base_size = 10) + ggtitle('F-half score')

acc; prc; roc; f1s; fhs
dev.off()


#################################
# make performance figures 
df_long = df %>% pivot_longer(cols = c('auPRC','f1_score','fpr'), 
                              values_to ='value',  names_to = 'variable')

pdf('figures/figure_mouse_ml_model-CNN_testSet_performance.pdf', 
    width = 5.5, height = 4)
ggplot(data = df_long, aes(x = celltype, y = value)) + 
  geom_boxplot(aes(fill = celltype)) + 
  geom_jitter(width = .25, pch = 21, aes(fill = celltype)) + 
  ylim(c(0,NA)) + 
  scale_fill_manual(values = colors, labels = celltypes) + 
  facet_grid(variable~lab + tissue, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = 11)
dev.off()

