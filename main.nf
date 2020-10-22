#! /usr/bin/env nextflow

SRA=Channel.of("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589")
chr_list = Channel.of(1..22,'MT','X','Y')

process getFASTQ {
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
    output:
    file 'Homo_sapiens.GRCh38.101.chr.gtf' into gtf_file

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz 
    """

}

process mapFASTQ {
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
        --outSAMunmapped None \
        --outFilterMismatchNmax 4 \
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
/*
process statAnalysis {

}
*/