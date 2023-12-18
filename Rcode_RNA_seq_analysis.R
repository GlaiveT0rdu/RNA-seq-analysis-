# Set the working directory
dir <- "D:/Work/Master 2/Internship 3/R_script/Rna_seq/data_preprocessing/" 
setwd(dir) 

# Load required libraries
library(readr)
library(DESeq2)
library(clusterProfiler)
library(AnnotationDbi)

# G1 Stg1
## Corresponds to G1 T1 Stg1
Br1 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br01_count.tsv", header=T, comment.char="#")
Br2 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br40_count.tsv", header=T, comment.char="#")
Br3 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br18_count.tsv", header=T, comment.char="#")

## Corresponds to G1 T2 Stg1
Br4 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br49_count.tsv", header=T, comment.char="#")
Br5 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br66_count.tsv", header=T, comment.char="#")
Br6 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br83_count.tsv", header=T, comment.char="#")

# G2 Stg1
## Corresponds to G2 T1 Stg1
Br1 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br05_count.tsv", header=T, comment.char="#")
Br2 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br27_count.tsv", header=T, comment.char="#")
Br3 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br48_count.tsv", header=T, comment.char="#")

## Corresponds to G2 T2 Stg1
Br4 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br53_count.tsv", header=T, comment.char="#")
Br5 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br70_count.tsv", header=T, comment.char="#")
Br6 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br96_count.tsv", header=T, comment.char="#")

# G3 Stg1
## Corresponds to G3 T1 Stg1
Br1 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br14_count.tsv", header=T, comment.char="#")
Br2 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br31_count.tsv", header=T, comment.char="#")
Br3 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br44_count.tsv", header=T, comment.char="#")

## Corresponds to G3 T2 Stg1
Br4 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br57_count.tsv", header=T, comment.char="#")
Br5 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br79_count.tsv", header=T, comment.char="#")
Br6 <- read.delim("D:/Work/Master 2/Internship 3/R_script/Rna_seq/Data/Br92_count.tsv", header=T, comment.char="#")


# Adjust dataframes
Br1 <- Br1[, c(1,7)]
Br2 <- Br2[, c(1,7)]
Br3 <- Br3[, c(1,7)]
Br4 <- Br4[, c(1,7)]
Br5 <- Br5[, c(1,7)]
Br6 <- Br6[, c(1,7)]

# Merge data
merged_data <- Reduce(function(x, y) merge(x, y, by="Geneid"), list(Br1, Br2, Br3, Br4, Br5, Br6))

# Rename the first column
colnames(merged_data)[1] <- "Gene"

# Count matrix
countData <- merged_data[, -1]
rownames(countData) <- merged_data$Gene

# Design matrix
coldata <- data.frame(row.names = colnames(countData), condition = c("Grp1", "Grp1","Grp1","Grp2", "Grp2", "Grp2"))

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ condition)

# Perform differential analysis
dds <- DESeq(dds)

# Get the results
res <- results(dds)

# Add base means
# res$baseMean <- baseMeans(dds)

# Extract normalized counts for each sample
norm_counts <- counts(dds, normalized=TRUE)

# Merge results and normalized counts
final_df <- cbind(as.data.frame(res), as.data.frame(norm_counts))

# Filter genes based on adjusted p-value and log2 fold change
filtered_genes <- subset(final_df, padj < 0.05 & abs(log2FoldChange) > 1)

# Export the final DataFrame to a CSV file
write.csv(filtered_G1, file = "filtered_results_G1_Stg1.csv")


###############################################################################################################################################################

# Set the working directory
dir <- "D:/Work/Master 2/Internship 3/R_script/Rna_seq/2-data_triaging/" 
setwd(dir) 

# Load required libraries
library(tidyverse)
library(ape)
library(VennDiagram)

# Read filtered results for each genotype
filtered_G1 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/2-data_triaging/filtered_results_G1_Stg1.csv")
filtered_G2 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/2-data_triaging/filtered_results_G2_Stg1.csv")
filtered_G3 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/2-data_triaging/filtered_results_G3_Stg1.csv")

# Rename columns
name <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "G1", "G2", "G3", "G4", "G5", "G6")
colnames(filtered_G1) <- name
colnames(filtered_G2) <- name
colnames(filtered_G3) <- name

# Merge dataframes
G1_G2 <- merge(filtered_G1, filtered_G2, by = "gene")
G1_G3 <- merge(filtered_G1, filtered_G3, by = "gene")
G2_G3 <- merge(filtered_G2, filtered_G3, by = "gene")
G1_G2_G3 <- merge(G1_G2, filtered_G3, by = "gene") 

# Extract gene sets
set1 <- filtered_G1$gene
set2 <- filtered_G2$gene
set3 <- filtered_G3$gene

genotype1 <- filtered_G1$gene
genotype2 <- filtered_G2$gene
genotype3 <- filtered_G3$gene
G123 <- G1_G2_G3$gene

# Define the output file path for the PDF
output_file <- "venn_diagram.pdf"

# Start writing to the PDF file
# pdf(output_file)

# Generate a Venn diagram
venn.plot <- venn.diagram(
  x = list(
    "Genotype 1" = genotype1,
    "Genotype 2" = genotype2,
    "Genotype 3" = genotype3
  ),
  category.names = c("Genotype 1", "Genotype 2", "Genotype 3"),
  output = "plot",
  filename = NULL,
  output.device = "Rdevice",
  imagetype = "pdf",
  category.colours = c("green", "blue", "orange"),
  fill = c("green", "blue", "orange")
)

# Display the Venn diagram
grid.draw(venn.plot)

# Close the PDF file
dev.off()

# Load CDS data
load("D:/Work/Master 2/Internship 3/R_script/Rna_seq/2-data_triaging/CDS_brassica.RData")

# Extract gene names from CDS data
CDS_Bt$extracted_gene <- sapply(strsplit(as.character(CDS_Bt$attributes), ";"), function(x) {
  gene_entry <- grep("gene=", x, value = TRUE)
  sub("gene=", "", gene_entry)
})

# Filter rows based on gene vector
filtered_genes_G1 <- CDS_Bt[CDS_Bt$extracted_gene %in% genotype1, ]
filtered_genes_G2 <- CDS_Bt[CDS_Bt$extracted_gene %in% genotype2, ]
filtered_genes_G3 <- CDS_Bt[CDS_Bt$extracted_gene %in% genotype3, ]
filtered_genes_G123 <- CDS_Bt[CDS_Bt$extracted_gene %in% G123, ]

# Create an ID variable for each gene ID
pattern <- "(?<=ID=cds-).+(?=;Parent)"
pattern2 <- "(?<=product=).+(?=;protein_id)"

tmp <- filtered_genes_G3
tmp$attributes <- as.character(tmp$attributes)

tmp$ID <- str_match(tmp$attributes, pattern)
tmp$product <- str_match(tmp$attributes, pattern2)
tmp$ID <- as.character(tmp$ID)
filtered_genes_G3 <- tmp

save(filtered_genes_G1, filtered_genes_G2, filtered_genes_G3, filtered_genes_G123, file = "D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/BeforeBlast.RData")


# Read protein sequence data
Ro18 <- read.fasta(file = "D:/Work/Master 2/Internship 3/R_script/Rna_seq/2-data_triaging/GCF_000309985.2_CAAS_Brap_v3.01_protein.faa", seqtype = "AA", as.string = TRUE)

for (r in 1:nrow(tmp)) {
  tmp$AA_sequence[r] <- try(get(tmp$ID[r], Ro18), silent = TRUE)
}

tmp <- tmp[!grepl("Error", tmp$AA_sequence), ]

tmp <- tmp[, c(-9:-12, -14)]

tmp <- tmp[, c(9, 10)]
all_info <- tmp

aa_seq <- AAStringSet(all_info$AA_sequence)
names(aa_seq) <- all_info$ID

writeXStringSet(aa_seq, "RNA_seq_G3.fasta")

###############################################################################################################################################################

# Set the working directory
dir <- "D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/" 
setwd(dir) 

# Load required libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(ape)
library(plyr)

# Read data from CSV files
orthoG1 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/orthoG1.csv", header=FALSE, sep=";")
orthoG2 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/orthoG2.csv", header=FALSE, sep=";")
orthoG3 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/orthoG3.csv", header=FALSE, sep=";")
orthoG1G2G3 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/orthoG1G2G3.csv", header=FALSE, sep=";")

trouveG1 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/trouveG1.csv", header=FALSE, sep=";")
trouveG2 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/trouveG2.csv", header=FALSE, sep=";")
trouveG3 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/trouveG3.csv", header=FALSE, sep=";")
trouveG1G2G3 <- read.csv("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/trouveG1G2G3.csv", header=FALSE, sep=";")

load("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/output/BeforeBlast.RData")

# Combine dataframes
G1 <- cbind(filtered_genes_G1, trouveG1, orthoG1)
G2 <- cbind(filtered_genes_G2, trouveG2, orthoG2)
G123 <- cbind(filtered_genes_G123, trouveG1G2G3, orthoG1G2G3)


# Calculate the maximum length among the three dataframes
max_length <- max(nrow(filtered_genes_G3), nrow(trouveG3), nrow(orthoG3))

# Adjust the length of filtered_genes_G3 if necessary
if (nrow(filtered_genes_G3) < max_length) {
  diff_length <- max_length - nrow(filtered_genes_G3)
  extension <- as.data.frame(matrix(NA, ncol=ncol(filtered_genes_G3), nrow=diff_length))
  colnames(extension) <- colnames(filtered_genes_G3)
  filtered_genes_G3 <- rbind(filtered_genes_G3, extension)
}

# Adjust the length of trouveG3 if necessary
if (nrow(trouveG3) < max_length) {
  diff_length <- max_length - nrow(trouveG3)
  extension <- as.data.frame(matrix(NA, ncol=ncol(trouveG3), nrow=diff_length))
  colnames(extension) <- colnames(trouveG3)
  trouveG3 <- rbind(trouveG3, extension)
}

# Adjust the length of orthoG3 if necessary
if (nrow(orthoG3) < max_length) {
  diff_length <- max_length - nrow(orthoG3)
  extension <- as.data.frame(matrix(NA, ncol=ncol(orthoG3), nrow=diff_length))
  colnames(extension) <- colnames(orthoG3)
  orthoG3 <- rbind(orthoG3, extension)
}

# Combine the adjusted dataframes by columns
G3 <- cbind(filtered_genes_G3, trouveG3, orthoG3)

# Display the combined dataframe
print(G3)

# Extract product information
pattern <- "(?<=product=).+(?=;protein_id)"
G1$product <- str_match(G1$attributes, pattern)
G2$product <- str_match(G2$attributes, pattern)
G3$product <- str_match(G3$attributes, pattern)
G123$product <- str_match(G123$attributes, pattern)

# Select specific columns
G1 <- G1[, c(1, 14, 12, 4, 5, 10, 21, 15, 18)]
G2 <- G2[, c(1, 14, 12, 4, 5, 10, 21, 15, 18)]
G3 <- G3[, c(1, 13, 12, 4, 5, 10, 14, 17, 18)]
G123 <- G123[, c(1, 14, 12, 4, 5, 10, 21, 15, 18)]

# Set column names
name <- c("shotgun_sequence", "Reference_Sequence", "Turnip_protein", "start", "end", "chromosome", "function", "arabido_gene", "protein_name_arabido")
colnames(G1) <- name
colnames(G2) <- name
colnames(G3) <- name
colnames(G123) <- name

# Save dataframes
save(G1, G2, G3, G123, file = "D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/Finish.Rdata" )

###############################################################################################################################################################

# Set the working directory
dir <- "D:/Work/Master 2/Internship 3/R_script/Rna_seq/4-Functional_analysis/" 
setwd(dir) 

# Load required libraries
library(clusterProfiler)
library(org.At.tair.db)

# Load data from the previous step
load("D:/Work/Master 2/Internship 3/R_script/Rna_seq/3-Blast/Finish.Rdata")

# Extract Arabidopsis gene names
G1g <- G1$arabido_gene
G1g <- sub("\\.[0-9]+$", "", G1g)
G2g <- G2$arabido_gene
G2g <- sub("\\.[0.9]+$", "", G2g)
G3g <- G3$arabido_gene
G3g <- sub("\\.[0.9]+$", "", G3g)
G123g <- G123$arabido_gene
G123g <- sub("\\.[0.9]+$", "", G123g)



# Obtain unique Arabidopsis gene names from G123
tG123 <- unique(G123g)

# Convert gene names to ENTREZID and SYMBOL using org.At.tair.db
converted_genes <- bitr(tG123, fromType = "TAIR", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.At.tair.db)

# Perform Gene Ontology (GO) enrichment analysis
egoG123 <- enrichGO(gene          = converted_genes$ENTREZID,
                    OrgDb         = org.At.tair.db,
                    ont           = "ALL",           # Choose from "BP", "CC", "MF", "ALL"
                    pAdjustMethod = "BH",            # Benjamini-Hochberg correction
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)            # Get readable gene names instead of Entrez IDs

# Create a PDF file for the GO plot
pdf("D:/Work/Master 2/Internship 3/R_script/Rna_seq/4-Functional_analysis/GO_G2.pdf", width=10, height=10)

# Generate a dot plot for the top 30 significant categories
dotplot(egoG123, showCategory = 30)

# Close the PDF file
dev.off()

