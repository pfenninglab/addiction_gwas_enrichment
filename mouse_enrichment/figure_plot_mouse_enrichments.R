
# Set up stuff
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyr)
library(ggplot2)
library(ggplot2)
library(RColorBrewer)

query1 = c( "AgeOfInitiation", "CigarettesPerDay",  "SmokingInitiation", "SmokingCessation" ,
            "DrinksPerWeek",'Cannabis', "RiskyBehavior")
query2 = c('ChronicPain','CocaineDep','OpioidDep','OCD','CoffeePerDay')
query3 = c('Schizophrenia','EduAttain','SleepDuration')
query4 = c('BMD','CAD','LBM')





###############################################
# Mouse tissue/celltype halLifted to hg38 enrichments single cell THS
dir = 'ldsc_enrichments/Mouse_enrichments'
fileNames = list.files(path = dir, pattern = '.cell_type_results.txt', full.names = T)
names(fileNames) = ss(ss(ss(fileNames,'/', 3),'_',5),'-')
mouse = do.call('rbind',lapply(fileNames, read.delim))
mouse$phenotype = ss(rownames(mouse),'\\.')
tmp = unique(mouse$phenotype)
mouse$phenotype = factor(mouse$phenotype, c(query1, query2, query3, query4))
mouse$bigGWAS = mouse$phenotype %in% query1
mouse$relatedGWAS = mouse$phenotype %in% query2
mouse$other1GWAS = mouse$phenotype %in% query3
mouse$other2GWAS = mouse$phenotype %in% query4

# celltypes and subtypes
mouse$lab = factor(ss(mouse$Name,'_',1),levels = c('Mo2015','Pfenning'))
mouse$Celltype = ss(mouse$Name,'_',2)
mouse$Celltype = gsub('pos','+',mouse$Celltype)
mouse$Celltype = gsub('neg','-',mouse$Celltype)
mouse$Celltype = gsub('Pv','PV',mouse$Celltype)
mouse$Tissue = ss(mouse$Name,'_',3)
mouse$Tissue[mouse$Tissue=='Ctx'] = 'CTX'
mouse$Tissue[mouse$Tissue=='Cpu'] = 'CPU'
mouse$Tissue[mouse$Tissue=='Nac'] = 'NAc'
mouse$Tissue = factor(mouse$Tissue, levels = c('CTX','CPU','NAc'))

class = c('Bulk','INTACT','Mo2015')
celltypeClass = rev(c('EXC+', 'D1+', 'D2+', 'VIP+', 'PV+', 'SST+', 'bulk'))
mycolors = rev(brewer.pal(name="Set1", n = 7)) # 3 Mo et al. INTACT, 2 MSN celltypes, 1 sst, 1 bulk 
mouse$Celltype = factor(mouse$Celltype, levels = celltypeClass)

# log p-values for plotting
mouse$FDR = p.adjust(mouse$Coefficient_P_value,'fdr')
mouse$Log_FDR = -log10(mouse$FDR)
mouse$Log_P_value = -log10(mouse$Coefficient_P_value)
mouse$FDR_signif = factor(mouse$FDR < 0.05)
mouse = mouse[complete.cases(mouse),]

pdf('figure_plots/figure_mouse_addiction_enrichment_plots_20200803.pdf' , h = 5, w = 12)
ggplot(data= subset(mouse,bigGWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
  facet_grid(lab + Tissue~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
  ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(legend.position = 'bottom')  + 
  theme(axis.text.y = element_text(face="bold", size=12),
        strip.text.x = element_text(face="bold", size=11),
        strip.text.y = element_text(face="bold", size=11))
dev.off()

pdf('figure_plots/figure_mouse_relatedGWAS_enrichment_plots_20200803.pdf' , h = 5, w = 11)
ggplot(data= subset(mouse,relatedGWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
  facet_grid(lab + Tissue~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
  ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(legend.position = 'bottom')  + 
  theme(axis.text.y = element_text(face="bold", size=12),
        strip.text.x = element_text(face="bold", size=11),
        strip.text.y = element_text(face="bold", size=11))
dev.off()

pdf('figure_plots/figure_mouse_other1GWAS_enrichment_plots_20200803.pdf' , h = 5, w = 6)
ggplot(data= subset(mouse,other1GWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
  facet_grid(lab + Tissue~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
  ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(legend.position = 'bottom')  + 
  theme(axis.text.y = element_text(face="bold", size=12),
        strip.text.x = element_text(face="bold", size=11),
        strip.text.y = element_text(face="bold", size=11))
dev.off()

pdf('figure_plots/figure_mouse_other2GWAS_enrichment_plots_20200803.pdf' , h = 5, w = 5)
ggplot(data= subset(mouse,other2GWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
  geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
  scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
  scale_color_manual(values=c(NA,'black'), guide = 'none')+
  geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
  facet_grid(lab + Tissue~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
  ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
  theme(legend.position = 'bottom')  + 
  theme(axis.text.y = element_text(face="bold", size=12),
        strip.text.x = element_text(face="bold", size=11),
        strip.text.y = element_text(face="bold", size=11))
dev.off()





