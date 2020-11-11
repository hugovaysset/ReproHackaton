#! /usr/bin/env nextflow

Channel.fromPath('resources/metadata.csv').into {samples_metadata_1; samples_metadata_2}
SRA = Channel.of("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589") //Channel containing all the SRA id of the samples of interest
//source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa
chr_list = Channel.of(1..22,'MT','X','Y') //Channel containing all the human chromosomes (including mitochondrial DNA)

process getFASTQ {
//This process permits to collect the genomic sequence of the samples of interest 
//It creates two compressed fastqc files for each samples of interest (one corresponding to the 5'3' sequencing and the other to the 3'5' sequencing)

    input:
    val SRAID from SRA

    output:
    //we triplicate the channel to use it in three different processes : mapFASTQ, fastqc and fastq_screen
    tuple val("${SRAID}"), file("${SRAID}_1.fastq.gz"), file("${SRAID}_2.fastq.gz") into fastq_files_to_map, fastq_files_to_qc, fastq_files_to_screen
    
    script:
    """
    fastq-dump --gzip --split-files ${SRAID}
    """
}

process getChrSeq {
//This process permits to collect the genomic sequence of each human chromosome
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
//This process permits to index the entire human genom
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
//This process permits to collect the latest available version of the human genome annotations
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
//This process permits to align the samples of interest with the entire human genome
//It creates BAM files

    input:
    path ref from genome_idx
    tuple val(SRAID), file("read1.fa.gz"), file("read2.fa.gz") from fastq_files_to_map

    output:
    file "${SRAID}.bam" into bamfiles
    file "*Log.final.out" into star_stats

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
        --outfileNamePrefix ${SRAID}
        --limitBAMsortRAM ${task.memory.toBytes()} > ${SRAID}.bam 
    """
}

process indexBAM {
//This process permits to index the BAM files generated by the previous process mapFASTQ
//The BAM files indexation will speed processes that use them

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
//This process permits to obtain, for each sample of interest, genes level of expression
//It counts the number of reads aligned on each gene of the human genome

    input:
    file "input.gtf" from gtf_file
    file bam from indexedBAM_1.flatten().filter(~/.*bam$/).collect() 
    // as the indexedBAM channel contains both bam and the corresponding index, filter to only pass bam file to featureCounts


    output:
    file "output_gene.counts" into counts
    file "output_gene.counts.summary" into counts_summary

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a input.gtf -o output_gene.counts ${bam}
    """
}

process countExon {
//This process permits to obtain, for each sample of interest, genes level of expression in case of alternative splicing
//It counts the number of reads aligned on each exon 

    input:
    file "input.gtf" from gtf_file
    file bam from indexedBAM_2.flatten().filter(~/.*bam$/).collect() 
    // as the indexedBAM channel contains both bam and the corresponding index, filter to only pass bam file to featureCounts


    output:
    file "output_exon.counts" into exoncounts
    file "output_exon.counts.summary" into exoncounts_summary

    script:
    """
    featureCounts -T ${task.cpus} -f -s 0 -a input.gtf -o output_exon.counts ${bam}
    """
}

/*
process gtf_to_gff {
// This process permits to transform a gtf file to a gff file
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
//This process permits to perfom statistic analysis of the samples of interest
//DESEQ2 allows to determine if there is a significant difference in gene expression between the wild and the mutate condition
//It uses the counting of the number of reads per gene done in the previous process countFeature
    publishDir "results/DE_genes", mode: 'symlink'

    input:
    file input from counts  // output of featureCounts
    file metadata from samples_metadata_1  // coldata

    output:
    file "gene_express_FC.csv"

    script:
    """
    statsAnalysis.R ${input} ${metadata} "gene_express_FC.csv"
    """
}

process statAnalysisSplicing {
//This process permits to perfom statistic analysis of the samples of interest
//DEXseq allows to determine if there is a significant difference in gene expression in case of alternative splicing between the wild and the mutate condition
//It uses the counting of the number of reads per exon done in the previous process countExon
    publishDir "results/DE_splicing", mode: 'symlink'

    input:
    file input from exoncounts  // output of featureCounts
    file metadata from samples_metadata_2  // coldata
    file gtf from gtf_file

    output:
    file "DEXSeq_subread.rds"
    file "DEXseq_results*"

    script:
    """
    statsAnalysisSplicing.R ${input} ${metadata} ${gtf} "DEXSeq_subread.rds" ${task.cpus}
    """
}


//Quality control on reads :


process fastqc {
//Fastqc allows to evalue the quality of the reads made by RNA-sequencing
    publishDir "results/fastqc_results", mode: 'symlink'

    input:
    tuple val(SRAID), file(read1), file(read2) from fastq_files_to_qc

    output:
    file "${SRAID}_[1-2]_fastq.*" into fastqc_results

    script:
    """
    fastqc -f fastq -q --threads ${task.cpus} ${read1} ${read2}
    """

}

process get_fastq_screen_genomes {
// Get bowtie2 genome index used by fastq_screen
    output:
    path "FastQ_Screen_Genomes" into FastQ_Screen_Genomes

    script:
    """
    fastq_screen --get_genomes
    """
}

process fastq_screen {
//Fastq_screen is used to compare the reads to a set of sequence databases in order to detect potential samples contaminations
    publishDir "results/fqscreen_results", mode:'symlink'

    input:
    tuple val(SRAID), file(read1), file(read2) from fastq_files_to_screen
    path "FastQ_Screen_Genomes" from FastQ_Screen_Genomes

    output:
    file("*_screen.txt") into fastq_screen_txt
    file("*_screen.html")
    
    script:
    """
    sed -i "s|/.*/FastQ_Screen_Genomes|$PWD/FastQ_Screen_Genomes|g" ./FastQ_Screen_Genomes/fastq_screen.conf
    fastq_screen --threads ${task.cpus}\
                 --conf ./FastQ_Screen_Genomes/fastq_screen.conf
                 --aligner bowtie2\
                 ${read1} ${read2}
    """
}


process multi_qc{
//MultiQC summarises results from quality control modules (fastQC and FastQ Screen), mapping module (STAR) and featureCounts
//Return html file containing a graphical representation of these results
    publishDir "results/multiqc_results", mode:'symlink'

    input:
    file fastqc from fastqc_results.collect()
    file star from star_stats.collect()
    file 'outputs_gene.count.summary' from counts_summary
    file 'outputs_exon.count.summary' from exoncounts_summary
    file fastq_screen from fastq_screen_txt.collect()

    output:
    file "multiqc_report.html"
    file "multiqc_data"

    script:
    """
    multiqc .
    """
}
