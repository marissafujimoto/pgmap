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

# Install conda and give write permissions to conda folder
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && chmod 777 -R /opt/conda/

RUN mkdir -p /home/biodocker/bin
RUN mkdir /data /config

ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/home/biodocker/bin
ENV HOME=/home/biodocker

RUN conda config --add channels r
RUN conda config --add channels bioconda

RUN conda upgrade conda

################## Install fastxtool kit ##################

#RUN wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/#fastx_toolkit-0.0.14.tar.bz2
#RUN tar -xjvf fastx_toolkit-0.0.14.tar.bz2
#RUN cd fastx_toolkit-0.0.14 && ./configure && make


################## Install Bowtie2 ##################

ENV ZIP=bowtie2-2.4.1-linux-x86_64.zip
ENV URL=https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1
ENV FOLDER=bowtie2-2.4.1-linux-x86_64
ENV DST=/home/biodocker/bin
ENV ULOCAL=/usr/local/bin

RUN wget $URL/$ZIP/download -O $DST/$ZIP && \
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

#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

################## Install fastp ##################

RUN     apt-get update \
&&      apt-get install -y wget build-essential zlib1g-dev \
&&      cd /tmp \
&&      wget https://github.com/OpenGene/fastp/archive/v0.20.0.tar.gz\
&&      tar xf v0.20.0.tar.gz \
&&      cd fastp-0.20.0 \
&&      make \
&&      make install \
&&      cd / \
&&      rm -rf /tmp/* \
&&      apt-get autoremove -y wget build-essential zlib1g-dev \
&&      rm -rf /var/lib/apt/lists/*

################# COPY over config files 
COPY config config/

WORKDIR /data

