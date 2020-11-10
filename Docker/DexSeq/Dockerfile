FROM rocker/r-ver:4.0.3


RUN apt-get update && apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl git libopenmpi-dev openmpi-bin openmpi-doc libxml2-dev wget zlib1g-dev libbz2-dev liblzma-dev \
  && wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 \
  && tar xvfj htslib-1.11.tar.bz2 && cd htslib-1.11 \
  && ./configure && make && make install \
  && R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}' \
   && R -e 'BiocManager::install("DEXSeq")'
