# Base image
FROM rocker/tidyverse:4.3.2
LABEL maintainer="cansav09@gmail.com"

WORKDIR /rocker-build/

# Install apt-getable packages to start
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog
RUN apt-get install -y --no-install-recommends \
    libxt6 \
    libpoppler-cpp-dev \
    vim \
    libglpk40 \
    curl \
    gpg \
    bzip2

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    python3-pip \
    python3-dev 

## Put more python packages we need here 
################### Install Seqtk ###################
ENV ZIP=seqtk-1.4.tar.gz
ENV URL=https://github.com/lh3/seqtk/archive/refs/tags/v1.4.tar.gz
ENV FOLDER=seqtk-1.4
ENV DST=/home/seqtk

RUN mkdir -p $DST && \
    wget $URL -O $DST/$ZIP && \
    tar -xzvf $DST/$ZIP -C $DST && \
    rm $DST/$ZIP

WORKDIR $DST/$FOLDER

RUN make

################## Install Bowtie2 ##################

ENV ZIP=bowtie2-2.4.1-linux-x86_64.zip
ENV URL=https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1
ENV FOLDER=bowtie2-2.4.1-linux-x86_64
ENV DST=/home/biodocker/bin
ENV ULOCAL=/usr/local/bin

RUN mkdir -p $DST && \
    wget $URL/$ZIP/download -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

################## Install Samtools ##################

RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y build-essential wget \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/src

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make


################# COPY over config files 
COPY config /home/config/

#################### COPY scripts ###################
COPY scripts /home/scripts/

# TODO remove this once paths are not hard coded
#################### COPY data ######################
COPY data /home/data/

WORKDIR /data

