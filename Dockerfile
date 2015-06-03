FROM ubuntu:14.04
MAINTAINER Yu Peng <loneknightpy@gmail.com>
RUN apt-get update && apt-get install -y gcc autoconf automake g++ make
ADD . /root/idba
RUN cd /root/idba && ./build.sh
ENV PATH=${PATH}:/root/idba/bin/
