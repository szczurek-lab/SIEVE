FROM condaforge/mambaforge:latest as conda_bin

FROM mambaorg/micromamba:jammy

LABEL maintainer="Senbai Kang <kang@mimuw.edu.pl>"

ENV ENV_NAME sieve
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND noninteractive
ENV _JAVA_OPTIONS -Xmx10g

ARG JAVA_VERSION=8
ARG BEAST_2_VERSION_MAJOR=2.6
ARG BEAST_2_VERSION_MINOR=7
ARG SIEVE_VERSION=v0.15.6

USER root

# Copy conda binaries from condaforge/mambaforge:latest as it is required by snakemake --use-conda
COPY --from=conda_bin /opt/conda /opt/conda

WORKDIR /tmp

# Install system dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    software-properties-common \
    dirmngr \
    apt-utils \
    locales \
    wget \
    curl \
    vim \
    htop \
    unzip \
    ca-certificates \
    apt-transport-https \
    build-essential \
    autotools-dev \
    libicu-dev \
    libbz2-dev \
    libomp-dev \
    gsfonts \
    gnupg2 && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# Add zulu Java repository
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 0xB1998361219BD9C9 && \
    curl -O https://cdn.azul.com/zulu/bin/zulu-repo_1.0.0-3_all.deb && \
    dpkg -i zulu-repo_1.0.0-3_all.deb && \
    rm zulu-repo_1.0.0-3_all.deb

# Install system dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    git \
    zulu${JAVA_VERSION}-ca-jre && \
    apt clean

# Create a micromamba environment containing snakemake
COPY environment.yml environment.yml
RUN micromamba create -y -f environment.yml && \
    micromamba clean -afy && \
    rm -f environment.yml

WORKDIR /root/pkgs

# Install BEAST 2 and SIEVE
RUN wget -O BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz https://github.com/CompEvol/beast2/releases/download/v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}/BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz && \
    tar -xzf BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz && \
    rm -f BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz && \
    ln -s "${PWD}/beast/bin/beast" /usr/local/bin/beast && \
    ln -s "${PWD}/beast/bin/applauncher" /usr/local/bin/applauncher && \
    ln -s "${PWD}/beast/bin/packagemanager" /usr/local/bin/packagemanager && \
    packagemanager -add ORC && \
    mkdir -p /root/.beast/${BEAST_2_VERSION_MAJOR}/SIEVE && \
    cd /root/.beast/${BEAST_2_VERSION_MAJOR}/SIEVE && \
    wget -O SIEVE-${SIEVE_VERSION}.tar.gz https://codeload.github.com/szczurek-lab/SIEVE/tar.gz/refs/tags/${SIEVE_VERSION} && \
    tar -xzf SIEVE-${SIEVE_VERSION}.tar.gz && \
    rm -f SIEVE-${SIEVE_VERSION}.tar.gz && \
    cp SIEVE*/dist/SIEVE.${SIEVE_VERSION}.zip ${SIEVE_VERSION}.zip && \
    rm -rf SIEVE* && \
    unzip ${SIEVE_VERSION}.zip && \
    rm -f ${SIEVE_VERSION}.zip

VOLUME /root/data

WORKDIR /root/data

# ENTRYPOINT is provided in the base image of micromamba to activate the conda environment and is inherited here.
# CMD is used to provide the default parameters to the ENTRYPOINT.
CMD ["snakemake", "--cores", "all", "--use-conda", "--rerun-incomplete", "--keep-going", "--printshellcmds", "--rerun-triggers", "mtime"]
