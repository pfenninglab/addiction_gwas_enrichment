options(stringsAsFactors = F)

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

load('../gwas_enrichment/rdas/fuma_addiction_snps.rda')

indUnique = which(!duplicated(snpDat$rsID))
snpDat = snpDat[indUnique, ]

################################################################
# for 167bp  extract the reference and alternate alleles of addiction snps
offset <- (167-1)/2 # make for 167 bp sequence

# get fasta sequences for effect allele
seqEffect <- xscat(getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos-1),
                snpDat$effect_allele,
                getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos+1,snpDat$pos+offset),
                sep='')
names(seqEffect) = paste0(snpDat$rsID,'_effect')
# get fasta sequences non effect allele
seqNonEffect <- xscat(getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos-1),
                snpDat$non_effect_allele,
                getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos+1,snpDat$pos+offset),
                sep='')
names(seqNonEffect) = paste0(snpDat$rsID,'_non-effect')

# write ref/alt fasta of snps
writeXStringSet(seqEffect, 'FASTA/addiction_snps_allele_effect_167.fa', format = 'fasta',width = 500)
writeXStringSet(seqNonEffect, 'FASTA/addiction_snps_allele_nonEffect_167.fa', format = 'fasta',width = 500)




################################################################
# for 1000bp extract the reference and alternate alleles of addiction snps
# get fasta sequences for effect allele
offset <- (1000-1)/2 # make for 1000 bp sequence
seqEffect <- xscat(getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos-1),
                snpDat$effect_allele,
                getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos+1,snpDat$pos+offset),
                sep='')
names(seqEffect) = snpDat$rsID
# get fasta sequences non effect allele
seqNonEffect <- xscat(getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos-1),
                snpDat$non_effect_allele,
                getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos+1,snpDat$pos+offset),
                sep='')
names(seqNonEffect) = snpDat$rsID

# write ref/alt fasta of snps
writeXStringSet(seqEffect, 'FASTA/addiction_snps_allele_effect_1000.fa', format = 'fasta',width = 1000)
writeXStringSet(seqNonEffect, 'FASTA/addiction_snps_allele_nonEffect_1000.fa', format = 'fasta',width = 1000)


################################################################
# for 500bp extract the reference and alternate alleles of addiction snps
# get fasta sequences for effect allele
peakLength = 501
offset <- (peakLength-1)/2 # make for 501 bp sequence
seqEffect <- xscat(getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos-1),
                snpDat$effect_allele,
                getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos+1,snpDat$pos+offset),
                sep='')
names(seqEffect) = snpDat$rsID
# get fasta sequences non effect allele
seqNonEffect <- xscat(getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos-offset,snpDat$pos-1),
                snpDat$non_effect_allele,
                getSeq(genome,paste0('chr',snpDat$chr),snpDat$pos+1,snpDat$pos+offset),
                sep='')
names(seqNonEffect) = snpDat$rsID

# write ref/alt fasta of snps
writeXStringSet(seqEffect, 'FASTA/addiction_snps_allele_effect_501.fa', format = 'fasta',width = peakLength)
writeXStringSet(seqNonEffect, 'FASTA/addiction_snps_allele_nonEffect_501.fa', format = 'fasta',width = peakLength)






