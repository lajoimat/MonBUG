library(data.table)
library(GenomicRanges)

# Load Van Allen Dataset (https://www.ncbi.nlm.nih.gov/pubmed/26359337)
variants = fread("MNV_and_liftover_Sept27/data/mutations.txt")

head(variants)
table(variants$Variant_Type)

# Get read counts
variants$t_alt_count = floor(variants$t_count*variants$i_tumor_f)

# Keep SNV only
snv = variants[variants$Variant_Type == "SNP",]

# Create sample specific chromosomes
snv$sample_chrom = paste(snv$patient, snv$Chromosome, sep = ":")
snv

# Create GRanges object from our SNV
gr.snv = GRanges(snv$sample_chrom, IRanges(snv$Start_position, snv$End_position))
gr.snv

# Add variants fields
mcols(gr.snv) = snv
gr.snv

# Sort by sample's genomic positions
gr.snv = sort(gr.snv)

# Find and merge adjacent SNVs
gr.merged = reduce(gr.snv, with.revmap = TRUE)
gr.merged

# Update variant type
gr.merged$NEW_VARTYPE = ifelse(elementNROWS(gr.merged$revmap) > 1, "MNV", "SNV")
table(gr.merged$NEW_VARTYPE)
table(gr.merged$NEW_VARTYPE) / length(gr.merged)


# Get variant information back with "relist" function

# Patients
gr.merged$patient = relist(snv$patient[unlist(gr.merged$revmap)] , skeleton = gr.merged$revmap)
gr.merged[gr.merged$NEW_VARTYPE == "MNV"]

# Start positions
gr.merged$original_positions = relist(gr.snv$Start_position[unlist(gr.merged$revmap)] , skeleton = gr.merged$revmap)
gr.merged[gr.merged$NEW_VARTYPE == "MNV"]
gr.merged$original_positions = NULL

# Alleles
gr.merged$ref = relist( gr.snv$Reference_Allele[ unlist(gr.merged$revmap)], skeleton = gr.merged$revmap)
gr.merged$alt = relist( gr.snv$Tumor_Seq_Allele2[ unlist(gr.merged$revmap)], skeleton = gr.merged$revmap)
gr.merged[gr.merged$NEW_VARTYPE == "MNV"]

### Operations on CharacterList ###

# Remove duplicated patient id and convert back to Character
gr.merged$patient = as.character(unique(gr.merged$patient))
gr.merged[gr.merged$NEW_VARTYPE == "MNV"]

# Concatenate alleles
gr.merged$ref = as.character(unstrsplit(gr.merged$ref))
gr.merged$alt = as.character(unstrsplit(gr.merged$alt))

gr.merged[gr.merged$NEW_VARTYPE == "MNV"]



# Get min and max read counts
gr.merged$alt_count = relist( gr.snv$t_alt_count[ unlist(gr.merged$revmap)], skeleton = gr.merged$revmap)
gr.merged[gr.merged$NEW_VARTYPE == "MNV"]

min_alt_counts = min(gr.merged$alt_count[gr.merged$NEW_VARTYPE == "MNV"])
max_alt_counts = max(gr.merged$alt_count[gr.merged$NEW_VARTYPE == "MNV"])

# Check if everything makes sense
plot(min_alt_counts, max_alt_counts, main = "MNV read counts"); grid()
abline(a=0,b=1,col="red")



# Look at this BRAF V600M variant

# Gene
gr.merged$Hugo_Symbol = relist( gr.snv$Hugo_Symbol[ unlist(gr.merged$revmap)], skeleton = gr.merged$revmap)

# Protein_Change
gr.merged$Protein_Change = relist( gr.snv$Protein_Change[ unlist(gr.merged$revmap)], skeleton = gr.merged$revmap)

gr.merged[gr.merged$NEW_VARTYPE == "MNV" & sum(gr.merged$Hugo_Symbol %in% "BRAF") > 0]


