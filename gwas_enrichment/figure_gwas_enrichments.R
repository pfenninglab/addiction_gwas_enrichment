# Set up stuff
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

query1 = c( "AgeOfInitiation", "CigarettesPerDay",  "SmokingInitiation", "SmokingCessation" ,
            "DrinksPerWeek",'Cannabis', "RiskyBehavior")
query2 = c('ChronicPain','CocaineDep','OpioidDep','OCD','CoffeePerDay')
query3 = c('Schizophrenia','EduAttain','SleepDuration')
query4 = c('BMD','CAD','LBM')

################################################################################################
# Fullard et al. GWAS enrichments human open chromatin w/ Honeybadger background
regions = c("OFC", "VLPFC","DLPFC", "ACC", "INS", "STC","ITC","PMC", "PVC", "AMY","HIPP", "MDT", "NAC", "PUT")
dir = 'ldsc_enrichments/Fullard_enrichments/no_overlap'
fileNames = list.files(path = dir, pattern = '.cell_type_results.txt', full.names = T)
fileNames = fileNames[grepl('Fullard_Honeybadger',fileNames)]
names(fileNames) = ss(ss(basename(fileNames),'_',3),'-')

fullard = do.call('rbind',lapply(fileNames, read.delim))
fullard$phenotype = factor(ss(rownames(fullard),'\\.'))
tmp = levels(fullard$phenotype)
fullard$phenotype = factor(fullard$phenotype, c(query1, query2, query3, query4))
fullard$bigGWAS = fullard$phenotype %in% query1
fullard$relatedGWAS = fullard$phenotype %in% query2
fullard$other1GWAS = fullard$phenotype %in% query3
fullard$other2GWAS = fullard$phenotype %in% query4


# celltype
fullard$Region = factor(ss(fullard$Name,'_',1), rev(regions))
fullard$bigRegion = ifelse(fullard$Region %in% c("AMY","HIPP", "MDT", "NAC", "PUT"),
                           'Subcortex','Cortex')
fullard$Celltype = factor(ss(fullard$Name,'_',2))
fullard$Celltype = factor(fullard$Celltype,levels(fullard$Celltype))

# log p-values for plotting
fullard$FDR = p.adjust(fullard$Coefficient_P_value,'fdr')
fullard$P.bonferroni = p.adjust(fullard$Coefficient_P_value,'bonferroni')
fullard$Log_P_value = -log10(fullard$Coefficient_P_value)
fullard$Log_FDR = -log10(fullard$FDR)
fullard$FDR_signif = fullard$FDR < 0.05
fullard$Log_bonf = -log10(fullard$P.bonferroni)

ind = which(fullard$FDR_signif)
with(fullard[ind,], table(phenotype))
with(fullard[ind,], table(Region))

ind2 = which(fullard$P.bonferroni < 0.05)
with(fullard[ind2,], table(phenotype))
with(fullard[ind2,], table(Region))

# make pretty plots
pdf('figure_plots/figure_fullard_addiction_enrichment_plots_20210414.pdf' , h = 4.5, w = 11)
(fullard_addiction_gwas = ggplot(data= subset(fullard,bigGWAS==TRUE), 
                                 aes(x = Region, y = Log_FDR, fill = Celltype))  + 
    geom_point(aes(shape=Celltype, fill=Celltype, stroke = FDR_signif)) + 
    scale_shape_manual(values=c(21, 24))+
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
    facet_grid( bigRegion~phenotype, scales = 'free_y',space = 'free_y') +  coord_flip() + 
    ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')  + 
    theme(axis.text.y = element_text(face="bold", size=12),
          axis.text.x = element_text(size=9),
          strip.text.x = element_text(face="bold", size=10),
          strip.text.y = element_text(face="bold", size=11))) + 
    theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0))
    
dev.off()


pdf('figure_plots/figure_fullard_relatedGWAS_enrichment_plots_20200803.pdf' , h = 4.5, w = 11)
(fullard_addiction_gwas = ggplot(data= subset(fullard,relatedGWAS==TRUE), 
                                 aes(x = Region, y = Log_FDR, fill = Celltype))  + 
        geom_point(aes(shape=Celltype, fill=Celltype, stroke = FDR_signif)) + 
        scale_shape_manual(values=c(21, 24))+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid( bigRegion~phenotype, scales = 'free_y',space = 'free_y') +  coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'none')  + 
        theme(axis.text.y = element_text(face="bold", size=12), 
              strip.text.x = element_text(face="bold", size=12),
              strip.text.y = element_text(face="bold", size=12)))
dev.off()

pdf('figure_plots/figure_fullard_other1GWAS_enrichment_plots_20200803.pdf' , h = 6, w = 6)
(fullard_addiction_gwas = ggplot(data= subset(fullard,other1GWAS==TRUE), 
                                 aes(x = Region, y = Log_FDR, fill = Celltype))  + 
        geom_point(aes(shape=Celltype, fill=Celltype, stroke = FDR_signif)) + 
        scale_shape_manual(values=c(21, 24))+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid( bigRegion~phenotype, scales = 'free_y',space = 'free_y') +  coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'bottom')  + 
        theme(axis.text.y = element_text(face="bold", size=12), 
              strip.text.x = element_text(face="bold", size=12),
              strip.text.y = element_blank()))
dev.off()

pdf('figure_plots/figure_fullard_other2GWAS_enrichment_plots_20200803.pdf' , h = 6, w = 5)
(fullard_addiction_gwas = ggplot(data= subset(fullard,other2GWAS==TRUE), 
                                 aes(x = Region, y = Log_FDR, fill = Celltype))  + 
        geom_point(aes(shape=Celltype, fill=Celltype, stroke = FDR_signif)) + 
        scale_shape_manual(values=c(21, 24))+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid( bigRegion~phenotype, scales = 'free_y',space = 'free_y') +  coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'bottom')  + 
        theme(axis.text.y = element_text(face="bold", size=12), 
              strip.text.x = element_text(face="bold", size=12),
              strip.text.y = element_text(face="bold", size=12)))
dev.off()





###############################################
# Lake et al. GWAS enrichments single cell THS
dir = 'ldsc_enrichments/Lake_enrichments'
fileNames = list.files(path = dir, pattern = '.cell_type_results.txt', full.names = T)
names(fileNames) = ss(ss(ss(fileNames,'/', 3),'_',3),'-')
lake = do.call('rbind',lapply(fileNames, read.delim))
lake$phenotype = ss(rownames(lake),'\\.')
tmp = unique(lake$phenotype)
lake$phenotype = factor(lake$phenotype, c(query1, query2, query3, query4))
lake$bigGWAS = lake$phenotype %in% query1
lake$relatedGWAS = lake$phenotype %in% query2
lake$other1GWAS = lake$phenotype %in% query3
lake$other2GWAS = lake$phenotype %in% query4

# celltypes and subtypes
class = c('Celltype','Cell Subtype')
mainClass = c( 'Ex','In', 'Ast','End','Mic','Oli','Opc' )
lab2class = c( 'Cortical_Excitatory','Cortical_Inhibitory',
               'Astrocyte','Endothelia','Microglia','Oligo','OPC' )
names(lab2class) = mainClass
subClass = c( 'ExL23','ExL4','ExL56','InA', 'InB')
#mycolors = c( brewer.pal(name="Accent", n = 7),brewer.pal(name="Set1", n = 5)) #main vs. sub
# mycolors = rev(c('#8dd3c7','#ffffb3','#e41a1c','#ff7f00', '#fb8072', '#fdb462','#b3de69')) #main colors

### set the Corces 2020 cell type clusters 
lab2cell = c("Cortical_Excitatory","Cortical_Inhibitory", 
             "Hippocampal_Excitatory","Striatal_Inhibitory",
             "Nigral_Neurons", "Unclassified_Neurons",
             "Astrocyte",  "Cortical_Astrocyte", "Nigral_Astrocyte", 
             "Striatal_Astrocyte", "Microglia", 
             "Oligo", "OPC", "Nigral_OPC", 'Endothelia')

# pre-define the cell type colors
mycolors = c(pal_aaas("default", alpha = 1)(6), pal_jco("default", alpha = .75)(9))
names(mycolors) = lab2cell

lake$Celltype = lab2class[lake$Name]
lake$Celltype = factor(lake$Celltype,rev(lab2class))
lake$cell_group = factor(ifelse(grepl('Astro|Micro|OPC|Oligo|Endo', lake$Celltype), 
                           'Glia','Neuron'), c('Neuron','Glia'))
lake$Celltype = droplevels(lake$Celltype)
lake = lake[complete.cases(lake),]

# log p-values for plotting
lake$FDR = p.adjust(lake$Coefficient_P_value,'fdr')
lake$P.bonferroni = p.adjust(lake$Coefficient_P_value,'bonferroni')
lake$Log_P_value = -log10(lake$Coefficient_P_value)
lake$Log_FDR = -log10(lake$FDR)
lake$Log_bonf = -log10(lake$P.bonferroni)
lake$FDR_signif = lake$FDR < 0.05

pdf('figure_plots/figure_lake_addiction_enrichment_plots_20210216.pdf' , h = 3, w = 11)
(lake_addiction_gwas = ggplot(data= subset(lake,bigGWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
        geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
        scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
        scale_color_manual(values=c(NA,'black'), guide = 'none')+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'bottom')  + 
        theme(axis.text.y = element_text(size=10),
              strip.text.x = element_text(face="bold", size=9)))
dev.off()


pdf('figure_plots/figure_lake_relatedGWAS_enrichment_plots_20210216.pdf' , h = 3, w = 11)
(lake_addiction_gwas = ggplot(data= subset(lake,relatedGWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
        geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
        scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
        scale_color_manual(values=c(NA,'black'), guide = 'none')+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'bottom')  + 
        theme(axis.text.y = element_text( size=10),
              strip.text.x = element_text(face="bold", size=10)))
dev.off()



pdf('figure_plots/figure_lake_other1GWAS_enrichment_plots_20210216.pdf' , h = 3, w = 6)
(lake_addiction_gwas = ggplot(data= subset(lake,other1GWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
        geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
        scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
        scale_color_manual(values=c(NA,'black'), guide = 'none')+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'bottom')  + 
        theme(axis.text.y = element_text( size=10),
              strip.text.x = element_text(face="bold", size=10)))
dev.off()



pdf('figure_plots/figure_lake_other2GWAS_enrichment_plots_20210216.pdf' , h = 3, w = 5)
(lake_addiction_gwas = ggplot(data= subset(lake,other2GWAS==TRUE), aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
        geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
        scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
        scale_color_manual(values=c(NA,'black'), guide = 'none')+
        geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
        facet_grid(cell_group~phenotype, scales = 'free_y', space = 'free') + coord_flip() + 
        ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
        theme(legend.position = 'bottom')  + 
        theme(axis.text.y = element_text( size=10),
              strip.text.x = element_text(face="bold", size=10)))
dev.off()











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

pdf('figure_plots/figure_mouse_addiction_enrichment_plots_20200803.pdf' , h = 6, w = 11)
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
    geom_bar(stat = 'identity', position=position_dodge(), color = 'black') + 
    scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
    # scale_color_manual(values=c(NA,'black'), guide = 'none')+
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








###############################################
# Macaque tissue/celltype halLifted to hg38 enrichments single cell THS
dir = 'ldsc_enrichments/Macaque_scRNAseq_genemarker_enrichments/allgenes_background'
fileNames = list.files(path = dir, pattern = '.cell_type_results.txt', full.names = T)
names(fileNames) = ss(ss(basename(fileNames),'_',4),'-')
macaque = do.call('rbind',lapply(fileNames, read.delim))
macaque$phenotype = ss(rownames(macaque),'\\.')
tmp = unique(macaque$phenotype)
macaque$phenotype = factor(macaque$phenotype, c(query1, query2, query3, query4))
macaque$bigGWAS = macaque$phenotype %in% query1
macaque$relatedGWAS = macaque$phenotype %in% query2
macaque$other1GWAS = macaque$phenotype %in% query3
macaque$other2GWAS = macaque$phenotype %in% query4

macaque$Name = gsub('NUDAP','RXFP1',macaque$Name)
macaque$Name = gsub('Olfactory_Tubercle','CPNE4',macaque$Name)

# celltypes and subtypes
celltypes = c( 'D1_MSNs_common','D1_MSN_CPNE4','D1_MSN_RXFP1','D2_MSNs','Interneurons',
               'Astrocytes','Endothelial','Ependymal','Microglia', 'Mural','Oligodendrocytes','OPC')
# mycolors = rev(c( brewer.pal(name="Set1", n = 5),brewer.pal(name="Set3", n = 7))) #neuronal vs. non-neuronal
mycolors = rev(c('#377eb8','#d0d1e6','#74a9cf','#4daf4a','#ff7f00',brewer.pal(name="Set3", n = 7))) #neuronal vs. non-neuronal
macaque$Celltype = ss(macaque$Name,'_Cells')
# macaque$Celltype[macaque$Celltype=='D1_MSN_Olfactory_Tubercle'] = 'D1_MSN_OlfTub'
macaque$Celltype[macaque$Celltype=='Oligodendrocyte_Precursors'] = 'OPC'
macaque$Celltype =factor(macaque$Celltype,rev(celltypes))
macaque$Group = factor(ifelse(grepl('MSN|Inter',macaque$Celltype),'Neuronal','Glial'),c('Neuronal','Glial'))

# log p-values for plotting
macaque$FDR = p.adjust(macaque$Coefficient_P_value,'fdr')
macaque$Log_FDR = -log10(macaque$FDR)
macaque$FDR_signif = macaque$FDR < 0.05

pdf('figure_plots/figure_macaque_addiction_enrichment_plots_allgenes_20200613.pdf' , h = 5, w = 11)
ggplot(data= subset(macaque,bigGWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
    geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
    scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
    scale_color_manual(values=c(NA,'black'), guide = 'none')+
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
    facet_grid(Group~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
    ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')  + 
    theme(axis.text.y = element_text(face="bold", size=12),
          strip.text.x = element_text(face="bold", size=10),
          strip.text.y = element_text(face="bold", size=11))
dev.off()

pdf('figure_plots/figure_macaque_relatedGWAS_enrichment_plots_allgenes_20200613.pdf' , h = 4, w = 11)
ggplot(data= subset(macaque,relatedGWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
    geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
    scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
    scale_color_manual(values=c(NA,'black'), guide = 'none')+
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
    facet_grid(Group~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
    ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')  + 
    theme(axis.text.y = element_text(face="bold", size=12),
          strip.text.x = element_text(face="bold", size=11),
          strip.text.y = element_text(face="bold", size=11))
dev.off()

pdf('figure_plots/figure_macaque_other1GWAS_enrichment_plots_allgenes_20200613.pdf' , h = 4, w = 6)
ggplot(data= subset(macaque,other1GWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
    geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
    scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
    scale_color_manual(values=c(NA,'black'), guide = 'none')+
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
    facet_grid(Group~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
    ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')  + 
    theme(axis.text.y = element_text(face="bold", size=12),
          strip.text.x = element_text(face="bold", size=11),
          strip.text.y = element_text(face="bold", size=11))
dev.off()

pdf('figure_plots/figure_macaque_other2GWAS_enrichment_plots_allgenes_20200613.pdf' , h = 4, w = 5)
ggplot(data= subset(macaque,other2GWAS==TRUE), 
       aes(x = Celltype, y = Log_FDR, fill = Celltype))  + 
    geom_bar(stat = 'identity', position=position_dodge(),aes(color = FDR_signif)) + 
    scale_fill_manual(values = mycolors, name = "Celltype", guide = 'none')+
    scale_color_manual(values=c(NA,'black'), guide = 'none')+
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'dashed') + 
    facet_grid(Group~phenotype, scales = 'free_y',space = 'free_y') + coord_flip() + 
    ylab('-log10(FDR)') + theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')  + 
    theme(axis.text.y = element_text(face="bold", size=12),
          strip.text.x = element_text(face="bold", size=11),
          strip.text.y = element_text(face="bold", size=11))
dev.off()




