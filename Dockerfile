FROM ubuntu:15.10
MAINTAINER Yu Peng <loneknightpy@gmail.com>
RUN apt-get update
RUN apt-get install -y openjdk-8-jdk
RUN apt-get install -y pkg-config zip g++ zlib1g-dev unzip bash-completion wget
RUN apt-get install -y gcc autoconf automake g++ make
RUN wget https://github.com/bazelbuild/bazel/releases/download/0.2.0/bazel_0.2.0-linux-x86_64.deb
RUN dpkg -i bazel_0.2.0-linux-x86_64.deb
ADD . /root/idba
RUN apt-get install -y vim
RUN cd /root/idba && bazel build //:*
