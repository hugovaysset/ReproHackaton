FROM ubuntu:18.04

RUN apt-get update \
	&& apt-get install -y wget libxml-libxml-perl uuid-runtime

RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz \
    && tar -vxzf sratoolkit.tar.gz -C /opt/ && rm sratoolkit.tar.gz
RUN mkdir /SRADATA
RUN mkdir /root/.ncbi

ENV PATH="/opt/sratoolkit.2.10.8-ubuntu64/bin/:${PATH}"

RUN printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg
RUN vdb-config --restore-defaults
RUN echo '/repository/user/main/public/root = "/SRADATA"' >> /root/.ncbi/user-settings.mkfg