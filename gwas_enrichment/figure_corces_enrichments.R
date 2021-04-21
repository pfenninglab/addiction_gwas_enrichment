ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(wesanderson)

### the different GWAS
query1 = c( "AgeOfInitiation", "CigarettesPerDay",  "SmokingInitiation", 
"SmokingCessation", "DrinksPerWeek",'Cannabis', "RiskyBehavior")
query2 = c('ChronicPain','CocaineDep','OpioidDep','OCD','CoffeePerDay')
query3 = c('Schizophrenia','EduAttain','SleepDuration')
query4 = c('BMD','CAD','LBM')

### read in the cell type enrichments
enrichments_fn = list.files(path = 'ldsc_enrichments/Corces_enrichments',
    pattern = '.cell_type_results.txt', full.names = T)
names(enrichments_fn) = basename(enrichments_fn)
input = lapply(enrichments_fn, read_tsv) %>% bind_rows(.id = 'file')
input$Name %>% table()

### set the Corces 2020 cell type clusters 
lab2cell = c("Cortical_Excitatory","Cortical_Inhibitory", 
             "Hippocampal_Excitatory","Striatal_Inhibitory",
             "Nigral_Neurons", "Unclassified_Neurons",
             "Astrocyte",  "Cortical_Astrocyte", "Nigral_Astrocyte", 
             "Striatal_Astrocyte", "Microglia", 
             "Oligo", "OPC", "Nigral_OPC")
names(lab2cell) = c( 
             "Isocortical_excitatory", "Isocortical_inhibitory", 
             "Hippocampal_excitatory", "Striatal_inhibitory",
             "Nigral_neurons",  "Neurons",
             "Astrocyte",  "Isocortical_astrocyte", "Nigral_astrocyte", 
             "Striatal_astrocyte", "Microglia", 
             "Oligodendrocyte", "OPC","Nigral_OPC")

# pre-define the cell type colors
mycolors = c(pal_aaas("default", alpha = 1)(6), pal_jco("default", alpha = .75)(8))
names(mycolors) = lab2cell

### prepare the dataframe
alpha = 0.05
df = input %>% 
  mutate(
    Name = as.character(Name),
    file = gsub('Corces_single_cell_','', as.character(file)),
    phenotype = factor(ss(file,'-')),
    phenotype = factor(phenotype, c(query1, query2, query3, query4)),
    bigGWAS = phenotype %in% query1,
    relatedGWAS = phenotype %in% query2,
    other1GWAS = phenotype %in% query3,
    other2GWAS = phenotype %in% query4,
    Celltype = lab2cell[ss(Name, '_c', 1)],
    Celltype = factor(Celltype, rev(lab2cell)),
    cell_group = case_when(
      grepl('Astro|Micro|OPC|Oligo', Celltype) ~ 'Glia',
      TRUE ~ 'Neuron'
    ),
    cell_group = factor(cell_group, c('Neuron', 'Glia')),
    clusters = ss(Name, '_c', 2)
  ) %>%
  mutate(
    FDR = p.adjust(Coefficient_P_value, 'fdr'), 
    Log_FDR = -log10(FDR), 
    FDR_signif = FDR < alpha
  )

table(df$Celltype, df$cell_group)
table(df$Celltype, df$FDR_signif)

### Addiction traits
pdf('figure_plots/figure_corces_addiction_enrichment_plots_20210216.pdf' , h = 4, w = 11)
ggplot(data= subset(df,bigGWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(alpha), colour = 'red', linetype = 'dashed') + 
  facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + 
  coord_flip() + ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(axis.text.y = element_text(size=10),
        strip.text.x = element_text(face="bold", size=9))
dev.off()

### neuropsych related traits
pdf('figure_plots/figure_corces_relatedGWAS_enrichment_plots_20210216.pdf' , h = 4, w = 11)
ggplot(data= subset(df,relatedGWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(alpha), colour = 'red', linetype = 'dashed') + 
  facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + 
  coord_flip() + ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(legend.position = 'none')  + 
  theme(axis.text.y = element_text(size=10), 
        strip.text.x = element_text(face="bold", size=10),
        strip.text.y = element_text(face="bold", size=10))
dev.off()
 
pdf('figure_plots/figure_corces_other1GWAS_enrichment_plots_20210216.pdf', h = 4, w = 6)
ggplot(data= subset(df,other1GWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(alpha), colour = 'red', linetype = 'dashed') + 
  facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + 
  coord_flip() + ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(axis.text.y = element_text(size=10), 
        strip.text.x = element_text(face="bold", size=10),
        strip.text.y = element_blank())
dev.off()


pdf('figure_plots/figure_corces_other2GWAS_enrichment_plots_20210216.pdf', h = 4, w = 5)
ggplot(data= subset(df,other2GWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(alpha), colour = 'red', linetype = 'dashed') + 
  facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + 
  coord_flip() + ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(axis.text.y = element_text(size=10), 
        strip.text.x = element_text(face="bold", size=10),
        strip.text.y = element_text(face="bold", size=10))

dev.off()



