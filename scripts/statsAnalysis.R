#!/usr/bin/Rscript --slave

library("DESeq2")

# Get path arguments given in command line
args <- commandArgs(TRUE)

input_path <- args[1]
metadata_path <- args[2]  # generer les metadonnes dans un process a part ?
output_path <- args[3]

# Load data
countData <- as.matrix(read.csv(input_path, header=TRUE, sep="\t"))
metaData <- as.data.frame(read.csv(metadata_path, header=TRUE, sep=","))

# Instantiate the DESeq data class
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~ SF3B1, tidy=TRUE)
                              # SF3B1 can be WT or Mut : it's the tested factor

# Run analysis
dds <- DESeq(dds)

# Extract the results
res <- results(dds)
resOrdered <- res[order(res$padj),]
resFiltered <- subset(resOrdered, padj < 0.1)
write.csv( as.data.frame(resFiltered), file=output_path )

# Eventuellement travail de renommage et verification des colonnes
# we need : countdata (matrix of read counts) and coldata (metadata on the
# conditions, columns should be c("sample", "patient", "treatment", "time") )

# # Run analysis : 
# # The DESeq command wraps three steps
# # dds <- estimateSizeFactors(dds)
# # dds <- estimateDispersions(dds)
# # dds <- nbinomWaldTest(dds)