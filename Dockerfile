FROM ubuntu:20.04

RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev libboost-all-dev libeigen3-dev && \
	apt-get clean



WORKDIR /opt


### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
  bunzip2 htslib-1.17.tar.bz2 && \
  tar xvf htslib-1.17.tar && \
	cd htslib-1.17/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make -j 4 && \
	make install 
	
	 

RUN cp -r /usr/local/* /usr/
RUN wget https://github.com/ChaissonLab/aou_test/archive/refs/tags/test3.tar.gz
RUN tar zxvf test3.tar.gz
RUN cd aou_test-test3 && make CONDA_PREFIX=/usr/
RUN cd aou_test-test3 && cp ctyper2 /usr/bin

ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib"

