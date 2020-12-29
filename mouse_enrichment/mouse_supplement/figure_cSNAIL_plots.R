library(SummarizedExperiment)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(limma)

celltypeClass = rev(c('EXC', 'D1', 'D2', 'VIP', 'PV', 'SST', 'bulk'))
mycolors = rev(brewer.pal(name="Set1", n = 7)) # 3 Mo et al. INTACT, 2 MSN celltypes, 1 sst, 1 bulk 
names(mycolors) = celltypeClass

pd = readRDS('rdas/cSNAIL_phenotype_table_n31.RDS')
rse = readRDS('rdas/cSNAIL_featureCounts_RangeSummarizedExperiment_n31_20200827.RDS')
colData(rse)$Batch = pd$Batch
colData(rse)$Group = paste(colData(rse)$Tissue, colData(rse)$Celltype, sep = "_")
colData(rse)$Group = factor(colData(rse)$Group)

dds <- DESeqDataSet(rse, design = ~ Group)
vsd <- vst(dds, blind=FALSE)

pc = plotPCA(vsd, intgroup = c('Tissue','Celltype'), returnData = T)

pdf('plots/mouse_supplement_cSNAIL_pca.pdf', height = 4, width = 4.5)
ggplot(data = pc, aes(x = PC1, y = PC2, fill = Celltype, shape = Tissue)) + 
  geom_point(size = 3) + 
  scale_shape_manual( values = c(21, 22, 23))  + 
  scale_fill_manual(values= mycolors)  + 
  theme_bw(base_size  = 14) +
  guides(fill = guide_legend(override.aes=list(shape=21, size = 2.5)),
         shape = guide_legend(override.aes=list(size = 2.5))) + 
  theme(legend.text=element_text(size=10))
dev.off()



