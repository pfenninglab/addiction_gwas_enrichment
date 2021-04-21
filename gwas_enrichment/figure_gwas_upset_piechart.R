options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(ggplot2)
library(ggupset)
library(tidyverse, warn.conflicts = FALSE)
load(file = 'rdas/fuma_addiction_snps.rda')


pdf('plots/addiction_gwas_upsetplot_loci.pdf', height = 5, width = 5.5)
lociDat %>%
  as_tibble() %>% group_by(loci) %>%
  summarize(Phenotype = strsplit(phenotype, ';')) %>%
  ggplot(aes(x = Phenotype)) +
  geom_bar() + scale_x_upset(order_by = "degree") + 
  ggtitle('Shared & Unique Trait Loci') + 
  theme_bw(base_size = 12) + theme(plot.margin = margin(.5, .5, .5, 3.5, "cm"),
                                   axis.text.y = element_text(size=14))
dev.off()

pdf('plots/addiction_gwas_upsetplot_snp.pdf', height = 5, width = 5.5)
snpDat %>%
  as_tibble() %>% group_by(uniqID) %>%
  summarize(Phenotype = list(phenotype)) %>%
  ggplot(aes(x = Phenotype)) +
  geom_bar() + scale_x_upset(order_by = "degree") +
  #ggtitle('Pleiotropic Substance Use SNPs') +
  theme_bw(base_size = 12) + theme(plot.margin = margin(.5, .5, .5, 3.5, "cm"),
                                   axis.text.y = element_text(size=14))
dev.off()


#######################################
# read in GWAS snps annotated with FUMA
gwas  = c('AgeOfInitiation', 'CigarettesPerDay', 'SmokingInitiation', 'SmokingCessation',
          'DrinksPerWeek', 'Cannabis', 'RiskyBehavior')

fumaSnps = list.files(path = 'fuma_annotations', pattern = 'snps.txt', recursive = T, full.names = T)
names(fumaSnps) = gwas[sapply(gwas, grep, fumaSnps)]

#######################################
# keep snps with R^2 more than 0.8
snpList = lapply(fumaSnps, read.delim)
colNames = Reduce('intersect',lapply(snpList, names))
snpList = lapply(snpList, function(x) x = x[x$r2 >= 0.8,colNames] )
snpDat = do.call('rbind', snpList)
snpDat$phenotype = factor(ss(rownames(snpDat), '\\.'), levels = gwas)

#######################################
# plot SNP annotation type pie chart
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
func = c('upstream', 'downstream', 'upstream:downstream', 
         'intergenic', 'intronic', 
         'exonic', 'splicing', 'ncRNA','UTR5', 'UTR3',
         'unknown') # regroup ncRNA:XX into XX groups

mycolors = c(brewer.pal(5, 'Dark2'), brewer.pal(6, 'Pastel1'))

# regroup ncRNA:XX into XX groups
snpDat2 = snpDat
snpDat2$func[is.na(snpDat2$func)] = 'unknown'
snpDat2$func[grepl('exonic',snpDat2$func)] = 'exonic'
snpDat2$func[grepl('splicing',snpDat2$func)] = 'splicing'
snpDat2$func[grepl('intronic',snpDat2$func)] = 'intronic'
snpDat2$func[grepl('ncRNA',snpDat2$func)] = 'ncRNA'
snpDat2$func = factor(snpDat2$func, levels = func)

snpFunc = snpDat2 %>%
  group_by(phenotype, func) %>%
  summarise(count = n_distinct(rsID))

# calculate percentages of each group
snpFunc = snpFunc  %>% 
  group_by(phenotype) %>% 
  mutate(percent=count/sum(count)*100)

pdf('figure_plots/fuma_snps_functional_annotation_pie.pdf', 
    width = 11, height = 2.3)
# base pie chart gg object
gg = ggplot(snpFunc, aes(x = '', y = percent, fill = func)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") # make PIE

# format colors, labels
gg + scale_fill_manual(values = mycolors) + 
  facet_wrap(~phenotype, nrow = 1) + xlab('Percent') + ylab("") +
  theme_bw(base_size = 11) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.position = "left",
        legend.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0), "cm"),
        strip.text.x = element_text(face= 'bold')) + 
  guides(fill = guide_legend(override.aes = list(size = 0.5))) + 
  theme(legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7)) + 
  theme(legend.title = element_blank())
dev.off()








