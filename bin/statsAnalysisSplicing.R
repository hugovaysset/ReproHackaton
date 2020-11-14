#!/usr/local/bin/Rscript --slave

library(DEXSeq)
library(magrittr)
library(stringr)

# Get path arguments given in command line
args <- commandArgs(TRUE)

input_path <- args[1]
metadata_path <- args[2]
gtf_path <- args[3]
output_path <- args[4]
cpus <- as.integer(args[5])

BPPARAM <-  BiocParallel::MulticoreParam(cpus)

# Parse GTF file for DEXSeq and remove duplicate entries
aggregates <- read.delim(gtf_path, stringsAsFactors = FALSE,
                         header = FALSE, comment.char = "#") %>%
  magrittr::set_colnames(c("chr", "source", "class", "start",
                           "end", "ex", "strand", "ex2", "attr")) %>%
  dplyr::filter(class == "exon") %>% # retain only exon entry
  dplyr::mutate(attr = str_remove_all(attr, "\"|=|;"),
                gene_id = sub(".*gene_id\\s(\\S+).*", "\\1", attr),
                exon_id = sub(".*exon_id\\s(\\S+).*", "\\1", attr),
                gene_name = sub(".*gene_name\\s(\\S+).*", "\\1", attr)) %>%
  dplyr::distinct(gene_id, exon_id, .keep_all = TRUE)

gene_exon <- dplyr::mutate(aggregates, join = str_c(gene_id, start, end, sep = "-")) %>%
  dplyr::select(join, exon_id, gene_name)

transcripts <- gsub(".*transcript_id\\s(\\S+).*", "\\1",
                    aggregates$attr)

exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start,
                                                          end = aggregates$end), 
                    strand = aggregates$strand, gene_name = aggregates$gene_name, 
                    ensemblID = aggregates$gene_id)

names(exoninfo) <- str_c(aggregates$gene_id, aggregates$exon_id, sep = ":")

names(transcripts) <- names(exoninfo)

metaData <- as.data.frame(read.csv(metadata_path, row.names = 1, sep = ","))
metaData$SF3B1 <- factor(metaData$SF3B1)

# Parse count and remove duplicate entries
read.table(input_path, header = TRUE) %>%
  dplyr::arrange(Geneid, Start, End) %>%
  dplyr::distinct() %>%
  dplyr::mutate(join = str_c(Geneid, Start, End, sep = "-")) %>%
  dplyr::inner_join(gene_exon) %>%
  dplyr::mutate(name = str_c(gene_id, exon_id, sep = ":")) %>%
  magrittr::set_rownames(.$name) %>%
  dplyr::select(-(Geneid:Length), -(join:name)) -> dcounts

## get genes and exon names out
splitted <- strsplit(rownames(dcounts), ":")
exons <- sapply(splitted, "[[", 2)
genesrle <- sapply(splitted, "[[", 1)

if (!all(rownames(dcounts) %in% names(exoninfo))) {
  stop("Count files do not correspond to the flattened annotation file")
}
matching <- match(rownames(dcounts), names(exoninfo))
stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))

design <- ~sample + exon + SF3B1:exon
dxd <- DEXSeqDataSet(dcounts, metaData, design, exons,
                     genesrle, exoninfo[matching], transcripts[matching])

dxr <- estimateSizeFactors(dxd) %>%
          estimateDispersions(BPPARAM = BPPARAM) %>%
          testForDEU(BPPARAM = BPPARAM) %>%
          estimateExonFoldChanges(BPPARAM = BPPARAM,fitExpToVar = "SF3B1") %>%
          DEXSeqResults()

# Plot if differential splicing
results <- dxr[which(dxr$padj < 0.01), ]

if (length(results@listData[[1]]) != 0) {
  # get top 10 genes
  gene_interest <- as.data.frame(results) %>% 
    dplyr::top_n(10,log2fold_WT_Mut) %>% 
    dplyr::pull(groupID)
  
  for (gene in gene_interest) {
    png(filename = paste0("DEXseq_results_",gene,".png"))
    plotDEXSeq(dxr, gene, fitExpToVar = "SF3B1", legend = TRUE, cex.axis = 1.2, cex = 1.3, lwd = 2)
    dev.off()
  }
  # write table
  results %>% as.data.frame() %>% write.csv("DEXSeq_results.csv")
  
} else {
  png(filename = "DEXseq_results.png")
  plot(0,0)
  text(0, 0, "No differential splicing")
  dev.off()
}
saveRDS(dxr, output_path)
