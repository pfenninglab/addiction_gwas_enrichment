# cd gene_tracks
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.basic.annotation.gff3.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.basic.annotation.gtf.gz
# cd ..
# gunzip gene_tracks/gencode.v32.basic.annotation.gtf.gz

gtfToGenePred -genePredExt -geneNameAsName2 gene_tracks/gencode.v32.basic.annotation.gtf gene_tracks/gencode.v32.basic.geneID.genePred
awk 'BEGIN{FS="\t"; OFS="\t"} {$1=$12; print}' < gene_tracks/gencode.v32.basic.geneID.genePred > gene_tracks/gencode.v32.basic.genePred
genePredToBed gene_tracks/gencode.v32.basic.genePred gene_tracks/gencode.v32.basic.bed
bedtools sort -i gene_tracks/gencode.v32.basic.bed > gene_tracks/gencode.v32.basic.sorted.bed


