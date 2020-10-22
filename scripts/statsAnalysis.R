#!/usr/bin/Rscript --slave

# library("DESeq2")

args <- commandArgs(TRUE)

input_path <- args[1]
output_path <- args[2]

cat("Input/output paths:\n")
print(input_path)
print(output_path)

# countData <- as.matrix(read.csv(input_path, header=TRUE, sep=";"))

# # eventuellement travail de renommage et verification des colonnes
# # we need : countdata (matrix of read counts) and coldata (metadata on the
# # conditions)

# # Instantiate the DESeq data class
# dds <- DESeqDataSetFromMatrix(countData=countData, 
#                               colData=metaData, 
#                               design=~ WT + mutated, tidy=TRUE)

# # run analysis. The DESeq command wraps three steps
# # dds <- estimateSizeFactors(dds)
# # dds <- estimateDispersions(dds)
# # dds <- nbinomWaldTest(dds)
# dds <- DESeq(dds)

# # extract the results
# res <- results(dds)
# resOrdered <- res[order(res$pvalue),]
# resFiltered <- subset(resOrdered, padj < 0.1)
# write.csv( as.data.frame(resFiltered), file=output_path )
