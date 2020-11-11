FROM continuumio/miniconda3:4.8.2 AS build

RUN conda install -y -c conda-forge conda-pack \
    && conda create -y -n fastq-screen -c bioconda fastq-screen=0.14.0-0 \ 
    && conda pack -n fastq-screen -o /tmp/fastq-screen.tar \
    && mkdir venv && cd /venv && tar xf /tmp/fastq-screen.tar \
    && rm /tmp/fastq-screen.tar && /venv/bin/conda-unpack 

FROM ubuntu:18.04

COPY --from=build /venv /venv

RUN apt-get update && apt-get install -y --no-install-recommends wget \
    && rm -rf /var/lib/apt/lists/* 

ENV PATH="/venv/bin/:${PATH}"

SHELL ["/bin/bash", "-c"]

