library(rtracklayer)
library(data.table)
library(GenomicRanges)

# INSTALLING THE FOLLOWING PACKAGE TAKES TIME, I RECOMMEND NOT DOING IT DURING PRESENTATION
library(BSgenome.Hsapiens.UCSC.hg38) 

# Load Van Allen Dataset (https://www.ncbi.nlm.nih.gov/pubmed/26359337)
variants = fread("MNV_and_liftover_Sept27/data/mutations.txt")

# Load chain file (http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)
chain = import("MNV_and_liftover_Sept27/data/hg19ToHg38.over.chain")
chain


# Keep SNV only
snv = variants[variants$Variant_Type == "SNP",]


# Create GRanges object from our SNV
gr.snv = GRanges(snv$Chromosome, IRanges(snv$Start_position, snv$End_position))

gr.snv$original_ref = snv$Reference_Allele
gr.snv$original_alt = snv$Tumor_Seq_Allele2
strand(gr.snv) = "+" # All variants were reported with respect to reference strand (+)
names(gr.snv) = as.character(gr.snv)
gr.snv


# Change chromosome style to UCSC
seqnames(gr.snv)
seqlevelsStyle(gr.snv) <- "UCSC"
seqnames(gr.snv)


# Liftover
lgr.lifted = liftOver(gr.snv,chain = chain)

# Look for number of unambiguous mapping
table(elementNROWS(lgr.lifted))

# Unlist 
gr.lifted = unlist(lgr.lifted)
gr.lifted


# Look for strand change
table(strand(gr.lifted))


# Keep strand change info
gr.lifted$strand_change = strand(gr.lifted) == "-"
gr.lifted


# Keep original position (as a key identifier)
gr.lifted$hg19_position = names(gr.lifted)
names(gr.lifted) = NULL # Just to make the print output clearer


# Get new Reference Allele
strand(gr.lifted) = "*"
gr.lifted$new_ref = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr.lifted))
gr.lifted


# Show new reference allele at strand changes
gr.lifted[gr.lifted$strand_change]
