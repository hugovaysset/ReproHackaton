#! /usr/bin/env nextflow

SRA="SRR628582"
chr_list = Channel.of(1..22,'M','X','Y')

process getFASTQ {
    input:
    val SRAID from SRA

    output:
    tuple file("${SRAID}_1.fastq.gz"), file("${SRAID}_2.fastq.gz") into fastq_files
    
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
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr${chr}.fa.gz
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

}

process mapFASTQ {

}

process indexBAM {

}

process countFeature {

}