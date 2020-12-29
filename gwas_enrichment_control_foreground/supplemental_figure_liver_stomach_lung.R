# Set up stuff
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyr)
library(ggplot2)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(wesanderson)


query1 = c( "AgeOfInitiation", "CigarettesPerDay",  "SmokingInitiation", "SmokingCessation" ,
            "DrinksPerWeek",'Cannabis', "RiskyBehavior")
query2 = c('ChronicPain','CocaineDep','OpioidDep','OCD','CoffeePerDay')
query3 = c('Schizophrenia','EduAttain','SleepDuration')
query4 = c('BMD','CAD','LBM')

###############################################
# opioid tissue/celltype halLifted to hg38 enrichments single cell THS
dir = 'Controls_EUR'
fileNames = list.files(path = dir, pattern = '.cell_type_results.txt', 
                       full.names = T, recursive = T)
names(fileNames) = paste0(ss(ss(basename(fileNames),'_',2),'-'),'_',
                          ss(fileNames,'/',2))
opioid = do.call('rbind',lapply(fileNames, read.delim))
opioid$phenotype = ss(ss(rownames(opioid),'\\.'),'_')
opioid$species = factor(ss(ss(rownames(opioid),'\\.'),'_',2),c('human','mouse'))
tmp = unique(opioid$phenotype)
opioid$phenotype = factor(opioid$phenotype, c(query1, query2, query3, query4))
opioid$bigGWAS = opioid$phenotype %in% query1
opioid$relatedGWAS = opioid$phenotype %in% query2
opioid$other1GWAS = opioid$phenotype %in% query3
opioid$other2GWAS = opioid$phenotype %in% query4
opioid$group = ifelse(opioid$bigGWAS, 'Addiction traits',
               ifelse(opioid$relatedGWAS,'Addiction-related traits',
               ifelse(opioid$other1GWAS,'Brain traits','Non-brain traits')))

# up vs down for each brain region
celltypes = c('preadipocyte','adipocyte', 'stomach','kidney','liver', 'lung')
opioid$Celltype = ss(opioid$Name,'_Cells')
opioid$Celltype =factor(opioid$Celltype,rev(celltypes))

# log p-values for plotting
opioid$FDR = p.adjust(opioid$Coefficient_P_value,'fdr')
opioid$Log_FDR = -log10(opioid$FDR)
opioid$FDR_signif = opioid$FDR < 0.05

pal <- wes_palette("Zissou1", 100, type = "continuous")

pdf('supplemental_fig_control_foregrounds_human_mouse.pdf' , h = 4.5, w = 11)
ggplot(data= opioid, aes(x = phenotype, y = Celltype, fill = Coefficient, alpha = FDR_signif))  +
    geom_tile(color = 'black', aes(size = FDR_signif)) +
    scale_fill_gradientn(colours = pal, name = 'LDSC Coefficient') + 
    scale_alpha_manual(values=c(.25,1), guide = 'none')+
    scale_size_manual(values=c(.1,1), guide = 'none')+
    facet_grid(species ~ group, scales = 'free',space = 'free') + 
    theme_bw(base_size = 12) +
    theme(legend.position = 'bottom')  +
    theme(axis.text.x = element_text(face="bold", size=12, angle = 30, hjust = 1),
          axis.text.y = element_text(face="bold", size=12),
          strip.text.x = element_text(face="bold", size=11),
          strip.text.y = element_text(face="bold", size=11),
          legend.text=element_text(size=10),
          legend.key.width=unit(1,"cm")) +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank())
dev.off()

