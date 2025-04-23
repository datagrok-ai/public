FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && \
    apt install -y --no-install-recommends curl bzip2 ca-certificates && \
    apt clean && rm -rf /var/lib/apt/lists/*

ENV MAMBA_EXE=/usr/local/bin/micromamba
ARG MAMBA_VERSION=1.5.3
RUN curl -L "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}" -o /tmp/micromamba.tar.bz2 && \
    mkdir -p /tmp/micromamba && \
    tar -xvjf /tmp/micromamba.tar.bz2 -C /tmp/micromamba --strip-components=1 && \
    mv /tmp/micromamba/micromamba "$MAMBA_EXE" && \
    rm -rf /tmp/micromamba /tmp/micromamba.tar.bz2

WORKDIR /app

COPY environment.yml .

RUN micromamba install -n base -c conda-forge awscli -y && \
    micromamba remove -n base awscli -y && \
    micromamba clean -afy

RUN micromamba env create -f environment.yml && \
    micromamba clean -afy

COPY . .

SHELL ["micromamba", "run", "-n", "myenv", "/bin/bash", "-c"]

EXPOSE 5000

ENTRYPOINT ["micromamba", "run", "-n", "myenv", "python", "grok_compute.py"]