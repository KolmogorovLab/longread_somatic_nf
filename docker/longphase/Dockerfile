FROM ubuntu:20.04
LABEL maintainer="https://github.com/twolinin/longphase"
LABEL version="1.7.3"

RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
	apt-get clean

WORKDIR /opt/longphase
RUN git clone https://github.com/twolinin/longphase.git /opt/longphase && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase/.git

### samtools
# 1.15
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar xvf samtools-1.15.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.15.tar.bz2 && \
	cd samtools-1.15/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools

ENV PATH="${PATH}":${HOME}/bin:/opt/longphase

CMD ["longphase", "phase", "--help"]
