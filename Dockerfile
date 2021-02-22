FROM snakemake/snakemake:v5.32.0


#### requirements for tabix and R
RUN apt-get update -y && apt-get install -y \
    apt-utils \
    bzip2 \
    gcc \
    make \
    ncurses-dev \
    wget \
    zlib1g-dev \
    r-base

RUN R -e "install.packages('doParallel',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && R -e "install.packages('foreach',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && R -e "install.packages('iterators',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && R -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"


RUN apt-get update && apt-get install -y --no-install-recommends build-essential wget unzip bzip2  libnss-sss && rm -rf /var/lib/apt/lists* \
    && conda install -c bioconda bedtools=2.27.0


ENV HTSLIB_INSTALL_DIR=/opt/htslib
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2

WORKDIR /tmp/htslib-1.3.2
RUN ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

RUN ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

ENTRYPOINT ["/bin/sh"]
