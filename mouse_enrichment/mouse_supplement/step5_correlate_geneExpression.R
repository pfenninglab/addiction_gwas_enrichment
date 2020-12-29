options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(SingleCellExperiment)
library(DESeq2)
library(IHW)
library(SummarizedExperiment)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(gridExtra)
library(reshape2)

celltypeClass = c('EXC', 'D1', 'D2', 'VIP', 'PV', 'SST', 'bulk')
mycolors = brewer.pal(name="Set1", n = 7)
names(mycolors) = celltypeClass

renameClusters = c('FC_1' = 'Ctx_VIP',
                 'FC_2' = 'Ctx_MGE',
                 'FC_6' = 'Ctx_EXC',
                 'STR_10' = 'Str_D1',
                 'STR_11' = 'Str_D2',
                 'STR_14' = 'Str_PV',
                 'STR_15' = 'Str_SST')

###########################################################################
## read in the summed expression across Saunders et al. cell type clusters
metaCells = readRDS('rdas/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS')
metaAnnot = readRDS('rdas/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS')
metaAnnot$tissue_cluster = ss(metaAnnot$tissue_subcluster,'-')

indKeep = with(metaAnnot, which(tissue %in% c('FC','STR') &
                                 class %in% c('NEURON') & 
                        tissue_cluster %in% names(renameClusters)))

metaCells = metaCells[,indKeep]
metaAnnot = metaAnnot[indKeep,]
metaAnnot$group = renameClusters[metaAnnot$tissue_cluster]

groupInd = split(seq(nrow(metaAnnot)), metaAnnot$group)
metaCell2 = do.call('cbind',lapply(groupInd, function(ind){
  tmp = metaCells[,ind]
  if(length(ind) > 1){
  toReturn = rowSums(tmp)
  } else{
    toReturn = tmp
    }
  return(toReturn)
}))
metaAnnot2 = data.frame(Group = names(groupInd),
                        Tissue = ss(names(groupInd),'_',1),
                        Celltype = gsub('Ctx_|Cpu_','',names(groupInd)))

# keep frontal cortex and striatum
seRNA = SummarizedExperiment(assays = list(counts = metaCell2), 
                             colData = metaAnnot2, 
                             rowData = data.frame(Symbol = rownames(metaCell2)))
colData(seRNA)$Assay = 'RNA'
vsdRNA = vst(DESeq2::DESeqDataSet(seRNA, ~1))

pcaRNA = plotPCA(vsdRNA, intgroup = c('Tissue','Celltype'), returnData = T)

pdf('plots/saunders_DropViz_FC_STR_pca.pdf', heigh = 4, width = 4.5)
ggplot(data = pcaRNA, aes(x = PC1, y = PC2, fill = Celltype, shape = Tissue)) + 
  geom_point(size = 3) + 
  scale_shape_manual( values = c(21, 22, 23))  + 
  scale_fill_manual(values= mycolors)  + 
  theme_bw(base_size  = 14) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()

###########################################################
## read in the TSS profiles from grouped ATAC-seq bigWigs
library(org.Mm.eg.db)
tabFiles = list.files(path = 'tables', pattern = '.tab.gz', full.names = T)
names(tabFiles) = ss(basename(tabFiles), '_compute')
tabFiles = tabFiles[c(grep('Ctx',names(tabFiles), value = T),
                      grep('D',names(tabFiles), value = T),
                      names(tabFiles)[!grepl('Ctx|D',names(tabFiles))])]

cn = c('seqname','start','end','Gene_ID','score','strand',seq(400))
seList = lapply(tabFiles, function(tab){
  tmp = read_tsv(tab, col_names = cn, skip = 1)
  tmp$Symbol = AnnotationDbi::mapIds(org.Mm.eg.db, keys = ss(tmp$Gene_ID,'\\.'), 
                                     column = "SYMBOL", keytype = "ENSEMBLTRANS", multiVals = "first")
  tmp = tmp[!is.na(tmp$Symbol) & !duplicated(tmp$Symbol),]
  se = SummarizedExperiment(assays = matrix(rowSums(tmp[,-c(1:6,407)]), ncol = 1),
                            rowData = data.frame(Symbol = tmp$Symbol))
  rm(tmp)
  return(se)
})

names(seList) = ss(basename(tabFiles), '_compute')
seList = seList
seATAC = do.call('cbind', seList)
names(assays(seATAC)) = 'counts'
rownames(seATAC) = rowData(seATAC)$Symbol
rowData(seATAC) = data.frame(Symbol = rowData(seATAC)[,'Symbol'])
colnames(seATAC) = names(seList) 
colData(seATAC)$Group = names(seList)
colData(seATAC)$Tissue = ss(names(seList), '_', 1)
colData(seATAC)$Celltype = ss(names(seList), '_', 2)
colData(seATAC)$Assay = 'ATAC'
assays(seATAC)[[1]] = round(assays(seATAC)[[1]] )
vsdATAC = vst(DESeq2::DESeqDataSet(seATAC, ~1))


#######################################
## keep genes & TSS in both datasets 
genesKeep = sort(intersect(rowData(seATAC)$Symbol, rowData(seRNA)$Symbol))
seATAC = seATAC[genesKeep,]
seRNA = seRNA[genesKeep,]

groups = c("Ctx_EXC","Ctx_PV", "Ctx_VIP", 
           "Cpu_D1", "Cpu_D2", "Cpu_PV", "Cpu_SST")

pdf('plots/correlation_plots_ATAC_RNA.pdf')
pList = lapply(groups, function(g){
  indATAC = match(g, colnames(seATAC))
  indRNA = match(g, colnames(seRNA))
  tmp = data.frame(RNA = log10(assays(seRNA)[[1]][genesKeep,indRNA] + 1),
                   ATAC = log10(assays(seATAC)[[1]][genesKeep,indATAC]+ 1 )  )
  p = ggplot(data = tmp, aes(x = ATAC, y = RNA)) + 
    geom_hex() + 
    ggtitle(g) + 
    xlab('log10(TSS counts)') + 
    ylab('log10(Gene counts)') + 
    theme_bw(base_size = 10)  + 
    theme(legend.position = 'none')
  p = p + stat_cor(aes(label = ..r.label.., size = 6), method = "pearson",
                   label.y.npc=1, label.x.npc = 0) 
  p = p + stat_cor(aes(label = ..r.label.., size = 6), method = "spearman",
                   cor.coef.name = 'rho',
                   label.y.npc=.8, label.x.npc = 0)
  p = ggplotGrob(p)
  return(p)
})
dev.off()

lay  <- rbind(c(1:2), c(3:4), c(5:6), c(NA, 7))
pdf('plots/correlation_plots_ATAC_RNA.pdf', height = 8, width = 3.5)
grid.arrange(grobs = pList, ncol = 2, layout_matrix=lay)
dev.off()





#############################
## make correlation matrix 
correlateMatrices = function(mat1, mat2, ind){
  lab1List = colnames(mat1)
  lab2List = colnames(mat2)
  
  # iterate over every column of mat1
  outList = lapply(lab1List, function(lab1){
    # iterate over every column of mat2
    outList = lapply(lab2List, function(lab2){
      # make sure getting same genes from both matrices
      vec1 = mat1[ind, lab1]
      vec2 = mat2[ind, lab2]
      # return Pearson correlation
      p = cor.test(vec1, vec2, method = 'pearson', exact=F)
      s = cor.test(vec1, vec2, method = 'spearman', exact=F)
      out = data.frame(lab1 = lab1, lab2 = lab2,
                       cor.pearson = p$estimate,
                       p.pearson = p$p.value,
                       p.spearman = s$p.value,
                       cor.spearman = s$estimate
                       )
      return(out)
    })
    outDat = do.call('rbind',outList)
    return(outDat)
  })

  outDat = do.call('rbind',outList)
  outDat$bonf.pearson = p.adjust(outDat$p.pearson, method = 'bonferroni')
  outDat$fdr.pearson =  p.adjust(outDat$p.pearson, method = 'fdr')

  outDat$bonf.spearman = p.adjust(outDat$p.spearman, method = 'bonferroni')
  outDat$fdr.spearman =  p.adjust(outDat$p.spearman, method = 'fdr')
  return(outDat)
}

corMat = correlateMatrices(mat1 = log10(assays(seRNA)[[1]] + 1),
                           mat2 = log10(assays(seATAC)[[1]] + 1 ),
                           ind = genesKeep)


#############################
## make ComplexHeatmap
mat = acast(corMat,  lab1 ~ lab2, value.var = 'cor.pearson')
mat = mat[,c(grep('Ctx',colnames(mat), value = T),
             grep('D',colnames(mat), value = T),
             colnames(mat)[!grepl('Ctx|D',colnames(mat))])]

mycolors2 = c(mycolors, 'MGE' = "#FF7F00")
tissue_col = c('Ctx' = 'palegreen', 'Str' = 'slategray', 
               'Cpu' = 'slategray', 'NAcc' = 'orchid')
col_fun = colorRamp2(c( .55, .65, .75), c("white", 'pink', "darkred"))
annot_side = rowAnnotation(Tissue = ss(rownames(mat), '_',1),
                           Cell_type = ss(rownames(mat), '_',2),
                    col = list(Cell_type = mycolors2, 
                               Tissue = tissue_col),
                    show_annotation_name = F)
annot_top = HeatmapAnnotation(Tissue = ss(colnames(mat), '_',1),
                              Cell_type = ss(colnames(mat), '_',2),
                        col = list(Cell_type = mycolors2,
                                   Tissue = tissue_col),
                        annotation_name_side = 'left')


  
lgd = Legend(col_fun = col_fun, )
ht3 = Heatmap(mat, name = "Correlation", top_annotation = annot_top, 
              right_annotation = annot_side, col = col_fun,
              cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F, 
              column_title = "cSNAIL ATAC Samples",
              # reorder columns by cell type then tissue
              column_order = order(ss(colnames(mat), '_',2), ss(colnames(mat), '_',1)),
              row_order = order(ss(rownames(mat), '_',2), ss(rownames(mat), '_',1)),
              row_title = 'Drop-seq')

pdf('plots/heatmap_plot_ATAC_RNA.pdf', width = 8.25, height = 2.5)
ht3
dev.off()


#############################
## make correlation matrix 
seBoth = cbind(seATAC, seRNA)
table(seBoth$Group, seBoth$Assay)

dds = DESeq2::DESeqDataSet(seBoth, ~Assay + Group)
vsd = vst(dds, blind = FALSE)
pca = plotPCA(vsd, intgroup = c('Tissue','Celltype','Assay','Group'), returnData = T)

pdf('plots/combine_RNA_ATAC_group_PCA.pdf', heigh = 4, width = 4.5)
ggplot(data = pca, aes(x = PC1, y = PC2, fill = Celltype, shape = Tissue)) + 
  geom_point(size = 3) + 
  scale_shape_manual( values = c(21, 22, 23))  + 
  scale_fill_manual(values= mycolors)  + 
  theme_bw(base_size  = 14) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()





