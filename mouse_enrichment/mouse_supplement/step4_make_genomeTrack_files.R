options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(reshape2)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(rtracklayer)

#############################
# get marker gene boundaries
genes = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mcols(genes)$symbol = mapIds(org.Mm.eg.db, keys = mcols(genes)$gene_id, 
                             column = "SYMBOL", keytype = "ENTREZID", 
                             multiVals = "first")

markers = c('Rbfox3', 'Snap25','Gad1', 'Slc17a7', 'Drd1', 
            'Drd2', 'Adora2a', 'Pvalb', 'Sst', 'Vip')
markerGenes = genes[match(markers, genes$symbol)]
tmp = width(markerGenes)
strand(markerGenes) = "*"
start(markerGenes) = start(markerGenes) - .1 * tmp
end(markerGenes) = end(markerGenes) + .1 * tmp
bedFile = 'bed/cSNAIL_markerGene_regions.bed'
export(markerGenes, con = bedFile)


###########################################################
# filter GTF file just around snpBlock, for fast plotting
gtf = import('bed/gencode.vM25.annotation.gtf.gz')

# filter out transcript_type and transcript_support_level
# https://www.gencodegenes.org/pages/data_format.html
gtf = gtf[gtf$transcript_support_level %in% c(1, 2)]
gtf = gtf[gtf$transcript_type %in% c('protein_coding', 'retained_intron', 
                                     'lncRNA','polymorphic_pseudogene')]

all(markers %in% unique(mcols(gtf)$gene_name))
gtf_filtered = gtf[mcols(gtf)$gene_name %in% markers]

# # get loci to plot around markerGenes
# # round to nearest 10kb
# start(markerGenes) = floor(start(markerGenes)/5/10^4 - 1)*5*10^4+1 - 10^5 
# # round to nearest 10kb
# end(markerGenes) = ceiling(end(markerGenes)/5/10^4 + 1)*5*10^4 + 10^5
# oo = findOverlaps(query = gtf, subject = markerGenes)
# gtf_filtered = gtf[unique(queryHits(oo))]
export(gtf_filtered, con = 'bed/TxDb.Mmusculus.UCSC.mm10.knownGene.filtered.gtf')



