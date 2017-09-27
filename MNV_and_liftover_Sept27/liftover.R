library(rtracklayer)
library(data.table)
library(GenomicRanges)

# Load Van Allen Dataset (https://www.ncbi.nlm.nih.gov/pubmed/26359337)
variants = fread("mutations.txt")


# Load chain file (http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)
chain = import("hg19ToHg38.over.chain")
chain
chain@listData


# Keep SNV only
snv = variants[variants$Variant_Type == "SNP",]


# Create GRanges object from our SNV
gr.snv = GRanges(snv$Chromosome, IRanges(snv$Start_position, snv$End_position))

gr.snv$ref = snv$Reference_Allele
gr.snv$alt = snv$Tumor_Seq_Allele2
strand(gr.snv) = "+" # All variants were reported with respect to reference strand (+)
names(gr.snv) = as.character(gr.snv)
gr.snv


# Change chromosome style to UCSC
seqnames(gr.snv)
seqlevelsStyle(gr.snv) <- "UCSC"
seqnames(gr.snv)


# Liftover
lgr.lifted = liftOver(gr.snv,chain = chain)
table(elementNROWS(lo))

gr.lifted = unlist(lgr.lifted)
gr.lifted

gr.lifted$hg19_position = names(gr.lifted)
table(strand(gr.lifted))
