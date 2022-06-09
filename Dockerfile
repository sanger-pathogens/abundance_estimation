FROM  gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/pathogens-base:0.2

ARG   DEBIAN_FRONTEND=noninteractive
ARG   HTSLIB_VERSION="1.9"
ARG   SAMTOOLS_VERSION="1.9"
ARG   BOWTIE2_VERSION="2.3.5"
ARG   INSTRAIN_VERSION="1.5.4"

ARG   SAMTOOLS_DIR="/usr/local/bin/samtools/${SAMTOOLS_VERSION}"
ARG   HTSLIB_DIR="/usr/local/bin/htslib/${HTSLIB_VERSION}"

WORKDIR /tmp

# Install required libraries
RUN   apt-get update  -qq -y && \
      apt-get install -qq -y apt-utils build-essential git zip zlib1g-dev libbz2-dev liblzma-dev libboost-dev && \
      apt-get install -qq -y libcurl4-gnutls-dev libssl-dev libbz2-dev libtool intltool libc-dev && \
      apt-get install -qq -y autoconf automake libncurses5-dev libncursesw5-dev libgd-dev libxml2-dev && \
      apt-get install --no-install-recommends -y unzip wget

# Install Python3
RUN   apt-get install --no-install-recommends -y python3 python3-distutils python3-pip python3-dev && \
      update-alternatives --install /usr/bin/python python /usr/bin/python3 2 && \
      update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 2

RUN   python3 --version && python --version && pip3 --version && pip --version

############## HTSLIB ################

RUN   wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
      tar -C /tmp/ -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
      cd /tmp/htslib-${HTSLIB_VERSION} && \
      autoheader && \
      autoconf && \
      ./configure --prefix="${HTSLIB_DIR}/" && \
      make && \
      make install prefix="${HTSLIB_DIR}/" && \
      rm -r /tmp/*

############# SAMTOOLS ################

RUN   wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
      tar -C /tmp/ -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
      cd samtools-${SAMTOOLS_VERSION} && \
      ./configure --prefix="${SAMTOOLS_DIR}/" && \
      make && \
      make install prefix="${SAMTOOLS_DIR}/" && \
      rm -r /tmp/*

ENV PATH=${SAMTOOLS_DIR}/bin:${PATH}

############# BOWTIE2 ################
RUN wget -q -O bowtie2.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip/download; \
	unzip bowtie2.zip -d /opt/; \
	ln -s /opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/ /opt/bowtie2; \
	rm bowtie2.zip
ENV PATH $PATH:/opt/bowtie2

############# INSTRAIN ################
RUN pip install boto3 && pip install instrain==${INSTRAIN_VERSION}

############# METAWRAP CUSTOM INSTALL ################
SHELL ["/bin/bash", "-c"]

# Install Miniconda package manager.
RUN wget -q -P /tmp \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

# Install mamba
RUN conda install -y -c conda-forge mamba

WORKDIR /opt
# download metaWRAP
RUN git clone https://github.com/bxlab/metaWRAP.git
ENV PATH="/opt/metaWRAP/bin:$PATH"

# custom trim-galore version
RUN sed -i 's|trim-galore 0.5.0|trim-galore 0.6.7|g' /opt/metaWRAP/conda_pkg/meta.yaml

# custom read_qc script
COPY metawrap_custom_files/read_qc.sh /opt/metaWRAP/bin/metawrap-modules

# custom config
COPY metawrap_custom_files/config-metawrap /opt/metaWRAP/bin

RUN mamba create -y -n metawrap-env python=2.7 && \
    source /opt/conda/etc/profile.d/conda.sh && \
    source ~/.bashrc && \
    conda activate metawrap-env && \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels ursky && \
    mamba install --only-deps -c ursky metawrap-mg

############ TIDY UP ###################
RUN   apt-get remove -y build-essential autoconf automake && \
      apt-get clean && \
      rm -rf /var/lib/apt/lists/*
