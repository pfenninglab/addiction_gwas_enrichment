# Set up stuff
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyverse)
library(ggplot2)

# read in the fisher overlap enrichment of functional annotations
annovar_fn = list.files(path = 'fuma_annotations', pattern = 'annov.stats.tx',
                        recursive = T, full.names = T)
names(annovar_fn) = annovar_fn %>% 
  ss(., '/', 2) %>% # folder name
  ss(., '_', 3) # GWAS phenotype name

annovar_df = do.call('rbind', lapply(annovar_fn, read.delim))
annovar_df$phenotype = rownames(annovar_df) %>% ss(.,'\\.')

# use Bonferroni correction of Fisher's exact P
annovar_df$p_bonferroni = annovar_df$fisher.P * 11

# show the data
with(annovar_df, table(phenotype, fisher.P < 0.01))
with(annovar_df, table(phenotype, p_bonferroni < 0.01))

# write the table out to xlsx
library(xlsx)

annovar_xlsx = file.path('tables','fuma_annovar_fisher_enrichment_gwas.xlsx')
xlsx::write.xlsx(annovar_df, file = annovar_xlsx)
