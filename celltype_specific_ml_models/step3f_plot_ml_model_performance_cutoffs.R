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
celltype = c('EXC','D1','D2')
colors = c('#e41a1c','#377eb8','#4daf4a')
names(colors) = celltype
tissue = c("Ctx", "Cpu", "Nac")
lab = c("Mo2015", "Pfenning")
groups = c('Ctx-EXC','Cpu-D1','Cpu-D2', 'Nac-D1', 'Nac-D2')

####################
# read in the data
df = read.delim('tables/SingleTask_IDR_bestModels_CNN_performance.tsv') %>%
  mutate(model_name = basename(model_name),
         group = factor(paste(tissue, celltype, sep = '-'), groups))

prediction_fn = list.files(path = 'predictions', pattern = '.predictions.txt',
                           full.names = T, recursive = T)
names(prediction_fn) = basename(prediction_fn)
input = prediction_fn %>% lapply(read_tsv, col_types = cols()) %>%
  data.table::rbindlist(idcol = 'file') %>%
  mutate( model_name = gsub('validation_negatives.predictions.txt|validation_positives.predictions.txt', 'h5', file),
          peaktype = case_when(grepl('positives', file) ~ 'pos', TRUE ~ 'neg'), 
          peaktype = factor(peaktype, c('pos', 'neg'))) %>% 
  right_join(df, by = 'model_name')

####################
# compute percent positives at different cutoffs
pred = input %>% group_by(model_name, peaktype, group, tissue, celltype) %>%
  summarise(cutoffs = seq(0, 1, .05),
            num = n(),
            pos.num = sapply(cutoffs, function(cut) sum(y_pred_score >= cut)),
            pos.prop = pos.num/ num) %>%
  ungroup()

pdf('figures/mouse_ml_model-CNN_calibration_cutoffs_20210413.pdf', width = 5.5, height = 2)
ggplot(pred, aes(x = cutoffs, y = pos.prop, fill = peaktype)) +
  geom_smooth(aes(color = peaktype), alpha = .5, se=FALSE) + 
  geom_boxplot(aes(group=paste(cutoffs, peaktype), fill = peaktype), 
               notch=FALSE, outlier.shape=NA,alpha = .8) + 
  scale_fill_manual(values = c('pos' = 'blue', 'neg' = 'red'), name = 'Validation set group') + 
  scale_color_manual(values = c('pos' = '#0000FF50', 'neg' = '#ff000050'), guide = 'none') + 
  xlab('Prediction Probability Threshold') + 
  ylab('Proportion Positive') + 
  theme_bw(base_size = 8) + facet_grid( ~ group) + 
  theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0))
dev.off()

tmp = input %>% filter(peaktype =='pos') %>% 
  group_by(model_name, group, id) %>% summarise(y_pred_score = mean(y_pred_score)) %>%
  ungroup() 

model_calibration = lapply(split(tmp$y_pred_score, tmp$group), ecdf)
calib_out_fn = 'rdas/SingleTask_IDR_CNN_validation_pos_calibration_ecdf.rds'
saveRDS(model_calibration, file = calib_out_fn)

