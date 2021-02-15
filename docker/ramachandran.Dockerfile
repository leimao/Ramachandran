FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

# Install package dependencies
RUN apt-get update
RUN apt-get install -y \
        wget \
        git \
        curl \
        language-pack-en \
        locales \
        locales-all \
        python3 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        python3-wheel \
        vim
RUN apt-get clean

RUN cd /usr/local/bin && \
    ln -s /usr/bin/python3 python && \
    ln -s /usr/bin/pip3 pip

# System locale
# Important for UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

RUN pip install numpy==1.20.1 \
                matplotlib==3.3.4 \
                aiohttp==3.7.3 \
                requests==2.25.1 \
                tqdm==4.56.2 \
                scipy==1.6.0 \
                twine==3.3.0
