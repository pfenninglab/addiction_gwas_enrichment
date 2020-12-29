# Set up stuff
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyverse)
library(ggplot2)
library(ggplot2)
library(RColorBrewer)

query1 = c( "AgeOfInitiation", "CigarettesPerDay",  "SmokingInitiation", "SmokingCessation" ,
            "DrinksPerWeek",'Cannabis', "RiskyBehavior")
query2 = c('ChronicPain','CocaineDep','OpioidDep','OCD','CoffeePerDay')
query3 = c('Schizophrenia','EduAttain','SleepDuration')
query4 = c('BMD','CAD','LBM')

################################################################################################
# Fullard et al. GWAS enrichments human open chromatin w/ Honeybadger background
regions = c("OFC", "VLPFC","DLPFC", "ACC", "INS", "STC","ITC","PMC", "PVC", "AMY","HIPP", "MDT", "PUT", "NAC")
regions2 = c("OFC", "VLPFC","DLPFC", "ACC", "STC", "PUT", "NAC")

dir = '../gwas_enrichment/ldsc_enrichments/Fullard_enrichments'
fileNames = list.files(path = dir, pattern = '.cell_type_results.txt', full.names = T)
fileNames = fileNames[grepl('Fullard_Honeybadger',fileNames)]
names(fileNames) = ss(ss(ss(fileNames,'/', 5),'_',3),'-')

fullard1 = do.call('rbind',lapply(fileNames, read.delim))
fullard1$phenotype = factor(ss(rownames(fullard1),'\\.'))
tmp = levels(fullard1$phenotype)
fullard1$phenotype = factor(fullard1$phenotype, c(query1, query2, query3, query4))
fullard1$bigGWAS = fullard1$phenotype %in% query1
fullard1$relatedGWAS = fullard1$phenotype %in% query2
fullard1$other1GWAS = fullard1$phenotype %in% query3
fullard1$other2GWAS = fullard1$phenotype %in% query4

# celltype
fullard1$Region = factor(ss(fullard1$Name,'_',1), rev(regions))
fullard1 = fullard1[fullard1$Region %in% regions2,]
fullard1$bigRegion = ifelse(fullard1$Region %in% c("AMY","HIPP", "MDT", "NAC", "PUT"),
                           'Subcortex','Cortex')
fullard1$Celltype = factor(ss(fullard1$Name,'_',2))
fullard1$Celltype = factor(fullard1$Celltype,levels(fullard1$Celltype))
fullard1$cutOff = 0.5
fullard1$Refinement = 'FACS'

#####################################################
# model-prioritization of Fullard et al. IDR peaks
dir2 = 'ldsc_enrichments/Fullard_cnn_scored_enrichments2'
fileNames2 = list.files(path = dir2, pattern = '.cell_type_results.txt', full.names = T)
fileNames2 = fileNames2[grepl('Fullard_CNNscored_enrichment',fileNames2)]
names(fileNames2) = ss(ss(ss(fileNames2,'/', 3),'_',4),'-')

fullard2 = do.call('rbind',lapply(fileNames2, read.delim))
fullard2$phenotype = factor(ss(rownames(fullard2),'\\.'))
tmp = levels(fullard2$phenotype)
fullard2$phenotype = factor(fullard2$phenotype, c(query1, query2, query3, query4))
fullard2$bigGWAS = fullard2$phenotype %in% query1
fullard2$relatedGWAS = fullard2$phenotype %in% query2
fullard2$other1GWAS = fullard2$phenotype %in% query3
fullard2$other2GWAS = fullard2$phenotype %in% query4

# celltype
fullard2$Region = factor(ss(fullard2$Name,'_',1), rev(regions))
fullard2 = fullard2[fullard2$Region %in% regions2,]
fullard2$bigRegion = ifelse(fullard2$Region %in% c("AMY","HIPP", "MDT", "PUT", "NAC"),
                            'Subcortex','Cortex')
fullard2$Celltype = factor(ss(fullard2$Name,'_',2))
fullard2$cutOff = as.numeric(ss(fullard2$Name,'_',3))
fullard2$Refinement = 'CNN'

#####################################################
# combine the original vs. the prioritized enrichments
# celltypes = c('EXC','PV','SST','VIP','D1','D2', 'NeuN+')
# colors = c('#e41a1c','#ff7f00','#ffff33','#984ea3','#377eb8','#4daf4a', 'black')
celltypes = c('EXC','D1','D2', 'NeuN+')
colors = c('#e41a1c','#377eb8','#4daf4a', 'black')
names(colors) = celltypes
fullard = rbind(fullard1,fullard2)
# fullard = fullard[fullard$cutOff==.75 | fullard$Celltype=='NeuN+' ,]
fullard = fullard[fullard$Celltype %in% celltypes,]
fullard$Region = factor(fullard$Region, rev(regions2))
fullard$Celltype = factor(fullard$Celltype, celltypes)
fullard = fullard[complete.cases(fullard),]

# log p-values for plotting
fullard$FDR = p.adjust(fullard$Coefficient_P_value,'fdr')
fullard$P.bonferroni = p.adjust(fullard$Coefficient_P_value,'bonferroni')
fullard$Log_P_value = -log10(fullard$Coefficient_P_value)
fullard$Log_FDR = -log10(fullard$FDR)
fullard$FDR_signif = fullard$FDR < 0.05
fullard$Log_bonf = -log10(fullard$P.bonferroni)

ind = 
with(fullard[ind,], table(phenotype))
with(fullard[ind,], table(Region))
with(fullard[ind,], table(Celltype))

fullard[which(fullard$FDR_signif &fullard$phenotype=='RiskyBehavior' ),]

ind2 = which(fullard$P.bonferroni < 0.05)
with(fullard[ind2,], table(phenotype))
with(fullard[ind2,], table(Region))


for (cut in unique(fullard$cutOff)){
  fdrPlotName = paste0('plots/figure_cnn_model_prioritized_fullard_addiction_enrichment_fdrPlot_',cut,'_20200730.pdf')
  # make pretty plots
  tmp = fullard %>% filter(bigGWAS==TRUE, cutOff==cut | Celltype=='NeuN+')
  pdf( fdrPlotName, h = 5, w = 8)
  p1 = ggplot(data= tmp, aes(x = Region, y = Log_FDR, fill = Celltype))  + 
    geom_jitter(aes(shape=Celltype, fill=Celltype, stroke = FDR_signif), size =2.5,width = 0.2) + 
    scale_shape_manual(values=c(24, 22, 23, 17))+
    scale_fill_manual(values=colors)+
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
    facet_grid( bigRegion~phenotype, scales = 'free_y',space = 'free_y') +  coord_flip() + 
    ylab('-log10(FDR)') + theme_bw(base_size = 12) + 
    theme(legend.position = 'bottom')  + 
    theme(axis.text.y = element_text(face="bold", size=10), 
          strip.text.x = element_text(face="bold", size=7),
          strip.text.y = element_text(face="bold", size=9))
  print(p1)
  dev.off()
  
  
  effPlotName = paste0('plots/figure_cnn_model_prioritized_fullard_addiction_enrichment_effectSizes_',cut,'_20200730.pdf')
  pdf(effPlotName , h = 4, w = 8)
  p2 = ggplot(data= tmp, aes(x = Region, y = Coefficient*10^8, fill = Celltype, alpha = FDR_signif))  + 
    geom_bar(stat = 'identity',aes(fill=Celltype, color = FDR_signif), position="dodge") + 
    geom_errorbar(aes(ymin= (Coefficient-Coefficient_std_error)*10^8, 
                      ymax=(Coefficient+Coefficient_std_error)*10^8), 
                  width = .3, position=position_dodge(width=.9)) +
    scale_fill_manual(values=colors, guide = FALSE)+
    scale_alpha_manual(values=c(.3,1), guide = FALSE)+
    scale_color_manual(values=c(NA,'black'), guide = FALSE)+
    geom_hline(yintercept = 0, colour = 'black') + 
    facet_grid( bigRegion~phenotype, scales = 'free',space = 'free_y') +  coord_flip() + 
    ylab(expression(paste('LDSC Coefficient (Scaled up by ', 10^8,")"))) + theme_bw(base_size = 12) + 
    theme(axis.text.y = element_text(face="bold", size=10), 
          strip.text.x = element_text(face="bold", size=7),
          strip.text.y = element_text(face="bold", size=9))
  print(p2)
  dev.off()
}









