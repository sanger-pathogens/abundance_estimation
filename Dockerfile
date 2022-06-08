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

############ TIDY UP ###################
RUN   apt-get remove -y build-essential autoconf automake && \
      apt-get clean && \
      rm -rf /var/lib/apt/lists/*
