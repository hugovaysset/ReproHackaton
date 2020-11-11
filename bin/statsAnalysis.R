#!/usr/bin/Rscript --slave

library("DESeq2")

# Get path arguments given in command line
args <- commandArgs(TRUE)

input_path <- args[1]
metadata_path <- args[2]  # generer les metadonnes dans un process a part ?
output_path <- args[3]

# Load data
countData <- read.csv(input_path, header=TRUE, comment.char="#", sep="	", row.names="Geneid")
metaData <- as.data.frame(read.csv(metadata_path, row.names=1, sep=","))
metaData$SF3B1 <- factor(metaData$SF3B1)
countData <- countData[, rownames(metaData)]

# Instantiate the DESeq data class
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~ SF3B1)
                              # SF3B1 can be WT or Mut : it's the tested factor

# Run analysis
dds <- DESeq(dds)

# Extract the results
res <- results(dds)
resOrdered <- res[order(res$padj),]
resFiltered <- subset(resOrdered, padj < 0.1)
write.csv( as.data.frame(resFiltered), file=output_path )

# # Run analysis : 
# # The DESeq command wraps three steps
# # dds <- estimateSizeFactors(dds)
# # dds <- estimateDispersions(dds)
# # dds <- nbinomWaldTest(dds)
