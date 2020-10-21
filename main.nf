#! /usr/bin/env nextflow

SRA="SRR628582"
chr_list = Channel.of(1..22,'MT','X','Y')

process getFASTQ {
    input:
    val SRAID from SRA

    output:
    tuple val(${SRAID}), file("${SRAID}_1.fastq.gz"), file("${SRAID}_2.fastq.gz") into fastq_files
    
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
        --limitBAMsortRAM ${task.memory} > ${SRAID}.bam
    """
}

process indexBAM {
    input:
    file bam from bamfiles

    output:
    file "${bam}.bai" into bamfilesindex

    script:
    """
    samtools index ${bam}
    """
}

process countFeature {
    input:
    file "input.gtf" from gtf_file
    file bam from bamfilesindex

    output:

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a input.gtf -o output.counts ${bam}
    """
}

process statAnalysis {

}