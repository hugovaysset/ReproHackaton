#! /usr/bin/env nextflow

SRA = Channel.of("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589") //Channel containing all the SRA id of the samples of interest
//source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa
chr_list = Channel.of(1..22,'MT','X','Y') //Channel containing all the human chromosomes (including mitochondrial DNA)

process getFASTQ { //This process permit to collect the genomic sequence of the samples of interest 
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

process getChrSeq { //This process permit to collect the genomic sequence of each human chromosome
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

process makeGenomeIndex { //This process permit to index the entire human genom
//It creates a fasta files (by combining all the fasta files generated at the getChrSeq process) that is then proceded by STAR which generated a lot of different files
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

process getAnnotations { //This process permit to collect the latest available version of the human genome annotations
//It creates a unique gtf file 
    output:
    file 'Homo_sapiens.GRCh38.101.chr.gtf' into gtf_file

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz 
    """

}

process mapFASTQ { //This process permit to align the samples of interest with the entire human genome
//It creates BAM files
    input:
    path ref from genome_idx
    tuple val(SRAID), file("read1.fa.gz"), file("read2.fa.gz") from fastq_files

    output:
    file "${SRAID}.bam" into bamfiles

    script:
    """
    STAR --genomeDir ${ref} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --readFilesIn read1.fa.gz read2.fa.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMunmapped None \             //Unmapped region are not keeped
        --outFilterMismatchNmax 4 \         //Mismatches are limited to 4
        --outFilterMultimapNmax 10 \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM ${task.memory.toBytes()} > ${SRAID}.bam 
    """
}

process indexBAM {
    input:
    file bam from bamfiles

    output:
    tuple file("${bam}.bai"), file("${bam}") into indexedBAM

    script:
    """
    samtools index ${bam}
    """
}

process countFeature {
    input:
    file "input.gtf" from gtf_file
    file bam from indexedBAM.flatten().filter(~/.*bam$/).collect() 
    // as the indexedBAM channel contains both bam and the corresponding index, filter to only pass bam file to featureCounts


    output:
    file "output.counts" into counts

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a input.gtf -o output.counts ${bam}
    """
}


process statAnalysis {
    // all paths are /tmp/ : it is the mount point of $baseDire in the container
    input:
    file input from fileChannel  // output of featureCounts
    val metadata from "/tmp/resources/metadata.csv"  // coldata
    val output from "/tmp/resources/gene_express_FC.csv"  //output path

    // no output required (end of the pipeline)
    // output:
    // file "$baseDir/resources/output.txt" into statsResults

    script:
    """
    /tmp/statsAnalysis.R ${input} ${metadata} ${output}
    """
}
