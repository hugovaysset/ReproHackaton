#! /usr/bin/env nextflow

Channel.fromPath('resources/metadata.csv').into {samples_metadata_1; samples_metadata_2}
SRA = Channel.of("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589") //Channel containing all the SRA id of the samples of interest
//source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa
chr_list = Channel.of(1..22,'MT','X','Y') //Channel containing all the human chromosomes (including mitochondrial DNA)

process getFASTQ {
//This process permit to collect the genomic sequence of the samples of interest 
//It creates two compressed fastqc files for each samples of interest (one corresponding to the 5'3' sequencing and the other to the 3'5' sequencing)

    input:
    val SRAID from SRA

    output:
    tuple val("${SRAID}"), file("${SRAID}_1.fastq.gz"), file("${SRAID}_2.fastq.gz") into fastq_files
    
    script:
    """
    fastq-dump --gzip --split-files ${SRAID}
    """
}

process getChrSeq {
//This process permit to collect the genomic sequence of each human chromosome
//It creates a compressed fasta file for each human chromosome

    input:
    val chr from chr_list
    
    output:
    file "*.fa.gz" into chrfa
    
    script:
    """
    wget -O chr${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    """
}

process makeGenomeIndex {
//This process permit to index the entire human genom
//The entire genome indexation will speeded processes that uses it
//It creates a fasta file (by combining all the fasta files generated at the getChrSeq process) that is then proceded by STAR which generated a lot of different files

    input:
    file '*.fa.gz' from chrfa.collect()

    output:
    path ref into genome_idx

    script:
    """
    gunzip -c *.fa.gz > ref.fa 
    mkdir ref
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ref.fa
    """
}

process getAnnotations {
//This process permit to collect the latest available version of the human genome annotations
//It creates a unique gtf file

    output:
    file 'Homo_sapiens.GRCh38.101.chr.gtf' into gtf_file

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz 
    """

}

process mapFASTQ {
//This process permit to align the samples of interest with the entire human genome
//It creates BAM files

    input:
    path ref from genome_idx
    tuple val(SRAID), file("read1.fa.gz"), file("read2.fa.gz") from fastq_files

    output:
    file "${SRAID}.bam" into bamfiles

//Unmapped region are not keeped and the number of mismatches are limited to 4
//Outpout files are forced to be in BAM format
    script:
    """
    STAR --genomeDir ${ref} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --readFilesIn read1.fa.gz read2.fa.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMunmapped None \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM ${task.memory.toBytes()} > ${SRAID}.bam 
    """
}

process indexBAM {
//This process permit to index the BAM files generated by the previous process mapFASTQ
//The BAM files indexation will speeded processes that uses them

    input:
    file bam from bamfiles

    output:
    tuple file("${bam}.bai"), file("${bam}") into indexedBAM_1, indexedBAM_2

    script:
    """
    samtools index ${bam}
    """
}

process countFeature {
//This process permit to obtain, for each sample of interest, genes level of expression
//It counts the number of reads aligned on each gene of the human genome

    input:
    file "input.gtf" from gtf_file
    file bam from indexedBAM_1.flatten().filter(~/.*bam$/).collect() 
    // as the indexedBAM channel contains both bam and the corresponding index, filter to only pass bam file to featureCounts


    output:
    file "output.counts" into counts

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a input.gtf -o output.counts ${bam}
    """
}

process countExon {
//This process permit to obtain, for each sample of interest, genes level of expression in case of alternative splicing
//It counts the number of reads aligned on each exon 

    input:
    file "input.gtf" from gtf_file
    file bam from indexedBAM_2.flatten().filter(~/.*bam$/).collect() 
    // as the indexedBAM channel contains both bam and the corresponding index, filter to only pass bam file to featureCounts


    output:
    file "output.counts" into exoncounts

    script:
    """
    featureCounts -T ${task.cpus} -f -s 0 -a input.gtf -o output.counts ${bam}
    """
}

/*
process gtf_to_gff {
// This process permit to transform a gtf file to a gff file
    input :
    file 'Annotation_Homo_sapiens.chr.gtf' from gtf_file

    output:
    file "Annotation_Homo_sapiens.chr.gff" into gff_file

    script:
    """
    dexseq_prepare_annotation.py Annotation_Homo_sapiens.chr.gtf Annotation_Homo_sapiens.chr.gff
    """
}
*/

process statAnalysis {
//This process permit to perfom statistic analysis of the samples of interest
//DESEQ2 allows to determine if there is a significant difference in gene expression between the wild and the mutate condition
//It uses the counting of the number of reads per gene done in the previous process countFeature
    publishDir "results/DE_genes", mode: 'symlink'

    input:
    file input from counts  // output of featureCounts
    file metadata from samples_metadata_1  // coldata

    // no output required (end of the pipeline)
    output:
    file "gene_express_FC.csv"

    script:
    """
    statsAnalysis.R ${input} ${metadata} "gene_express_FC.csv"
    """
}

process statAnalysisSplicing {
//This process permit to perfom statistic analysis of the samples of interest
//DEXseq allows to determine if there is a significant difference in gene expression in case of alternative splicing between the wild and the mutate condition
//It uses the counting of the number of reads per exon done in the previous process countExon
    publishDir "results/DE_splicing", mode: 'symlink'

    input:
    file input from exoncounts  // output of featureCounts
    file metadata from samples_metadata_2  // coldata
    file gtf from gtf_file

    // no output required (end of the pipeline)
    output:
    file "DEXSeq_subread.rds"
    file "DEXseq_results*"

    script:
    """
    statsAnalysisSplicing.R ${input} ${metadata} ${gtf} "DEXSeq_subread.rds" ${task.cpus}
    """
}
