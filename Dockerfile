FROM ubuntu:20.04 as builder

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install wget -y && \
    apt-get install gcc -y && \
    apt-get install cpp -y && \
    apt-get install make -y && \
    apt-get install parallel -y && \
    apt-get install cpanminus -y && \
    apt-get install curl -y && \
    apt-get install git -y && \ 
    apt-get install build-essential -y && \ 
    apt-get install zlib1g-dev -y && \
    apt-get install libncurses5-dev -y && \ 
    apt-get install libncursesw5-dev -y && \
    apt-get install libexpat1-dev -y && \
    apt-get install libdb-dev -y && \
    apt-get install locales -y && \
    apt-get install libssl-dev -y && \
    rm -r /var/lib/apt/lists/*


ENV PATH /opt/conda/bin:$PATH

RUN wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O miniconda.sh && \ 
    /bin/bash miniconda.sh -b -p /opt/conda/ && \
    rm miniconda.sh

RUN conda update -n base -c defaults conda -y && \
    conda install \
    -c conda-forge \
    -c bioconda \
    coreutils==8.31 \
    blat==35 \
    cap3==10.2011 \
    samtools==0.1.19 \
    -y && \
    conda clean --all -y

ENV PATH /opt/conda/bin:$PATH

RUN locale-gen en_US.UTF-8
RUN update-locale
ENV LC_ALL en_US.UTF-8

RUN umask 002
RUN mkdir -p /usr/local/perlbrew /root
ENV HOME /usr/local
ENV PERLBREW_ROOT /usr/local/perlbrew
ENV PERLBREW_HOME /usr/local/.perlbrew

RUN curl -kL http://install.perlbrew.pl | bash
ENV PATH /usr/local/perlbrew/bin:$PATH
ENV PERLBREW_PATH /usr/local/perlbrew/bin
ENV PERL_VERSION 5.10.1

RUN perlbrew --notest install 5.10.1
ENV PATH ${PERLBREW_ROOT}/perls/perl-$PERL_VERSION/bin:$PATH
RUN perlbrew switch perl-5.10.1
RUN perlbrew install-cpanm
RUN perlbrew info

RUN which perl
RUN perl -v

RUN which perl
RUN perl -v

RUN cpanm --force -i \
    Algorithm::Combinatorics@0.26

RUN cpanm --force -i \
    Set::IntSpan@1.19

RUN cpanm --force -i \ 
    enum \
    Data::Compare@1.22 \
    DBI@1.626 \
    Test::Deep@1.128 && \
    chown -R root:root /usr/local/.cpanm

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/archive/0.1.17.tar.gz && \
    tar -xzf 0.1.17.tar.gz && \
    cd samtools-0.1.17 && \
    sed -i "s/CFLAGS=.*$/CFLAGS= -g -Wall -O2 -fPIC/" Makefile && \
    cat Makefile && \
    make lib 

RUN SAMTOOLS="/tmp/samtools-0.1.17" cpanm --force -i Bio::DB::Sam@1.35 && chown -R root:root /usr/local/.cpanm

FROM ubuntu:20.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y libncursesw5 && \
    rm -r /var/lib/apt/lists/*

COPY --from=builder /opt/conda/bin/gfClient /opt/conda/bin/gfClient
COPY --from=builder /opt/conda/bin/gfServer /opt/conda/bin/gfServer
COPY --from=builder /opt/conda/bin/blat /opt/conda/bin/blat

COPY --from=builder /opt/conda/bin/cap3 /opt/conda/bin/cap3
COPY --from=builder /opt/conda/bin/samtools /opt/conda/bin/samtools
COPY --from=builder /opt/conda/lib /opt/conda/lib
COPY --from=builder /usr/bin/parallel /usr/bin/parallel

COPY --from=builder /usr/local/perlbrew/perls/perl-5.10.1 /usr/local/perlbrew/perls/perl-5.10.1
ENV PATH /usr/local/perlbrew/perls/perl-5.10.1/bin:$PATH

COPY src/scripts /opt/cicero/src/bin
COPY src/perllib /opt/cicero/src/perllib
COPY dependencies/lib/perl/* /opt/cicero/src/perllib/
COPY configs /opt/cicero/configs

ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/cicero/src/bin:${PATH}
ENV PERL5LIB /opt/cicero/src/perllib:${PERL5LIB}

ENTRYPOINT ["/opt/cicero/src/bin/Cicero.sh"]
CMD ["-h"]
