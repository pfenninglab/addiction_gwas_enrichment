# Set up stuff
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(ggplot2)
library(ggcorrplot)
library(tidyverse)
library(data.table)
library(reshape2)

query1 = c( "AgeOfInitiation", "CigarettesPerDay",  "SmokingInitiation", "SmokingCessation" ,
            "DrinksPerWeek",'Cannabis', "RiskyBehavior")
query2 = c('ChronicPain','CocaineDep','OpioidDep','OCD','CoffeePerDay')
query3 = c('Schizophrenia','EduAttain','SleepDuration')
query4 = c('BMD','CAD','LBM')


#####################
# read in the GWASs
sumstatsFile = '/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/mouse_enrichment/gwas_list_sumstats.txt'
df = read.delim(sumstatsFile, header = F, col.names = c('sumstatsFile'))
df$pheno = ss(basename(df$sumstatsFile),'-',1)

ldscCorrFiles = file.path('tables',paste0('p1_',ss(basename(df$sumstatsFile),'\\.',1),'_correlation.log'))
dfCor = do.call('rbind',lapply(ldscCorrFiles, fread, skip = 'Summary of Genetic Correlation Results'))
dfCor = distinct(dfCor)

dfCor$p1 = factor(ss(basename(dfCor$p1),'-',1), c(query1, query2, query3, query4))
dfCor$p2 = factor(ss(basename(dfCor$p2),'-',1), c(query1, query2, query3, query4))

dfCor$FDR = p.adjust(dfCor$p)

corMat = acast(subset(dfCor, p1 %in% query1 & p2 %in% query1), p1~p2, value.var = 'rg')
pMat = acast(subset(dfCor, p1 %in% query1 & p2 %in% query1), p1~p2, value.var = 'p')
fMat = acast(subset(dfCor, p1 %in% query1 & p2 %in% query1), p1~p2, value.var = 'FDR')

pdf('plots/addiction_GWAS_correlation.pdf', heigh = 5.5, width = 5.5)
ggcorrplot(
  corMat,
  p.mat = fMat, # correlation FDR cutoff
  hc.order = TRUE,
  type = "lower",
  # insig = "blank", 
  outline.color = "black",
  lab = TRUE, 
  legend.title = expression('r'['g']), 
  show.diag = TRUE,
  lab_size = 4,
  sig.level = 0.05
) + theme(axis.text.x = element_text( size=12, face="bold"),
          axis.text.y = element_text( size=12, face="bold")
)
dev.off()


