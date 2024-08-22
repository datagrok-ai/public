FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

SHELL ["/bin/bash","-c"]

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV MAMBA_EXE=/bin/micromamba

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl bzip2 ca-certificates ; \
     apt-get clean

RUN set -ex ; \
    curl -L "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}" --output /tmp/micromamba.tar.bz2; \
    mkdir -p /tmp/micromamba ; \
    tar --bzip2 --extract -C /tmp/micromamba --strip-components=1 -f /tmp/micromamba.tar.bz2 ; \
    mv /tmp/micromamba/micromamba "${MAMBA_EXE}" ; \
    rm -rf /tmp/micromamba

ARG MAMBA_USER=mambauser
ARG MAMBA_USER_ID=57439
ARG MAMBA_USER_GID=57439
ENV MAMBA_USER=$MAMBA_USER
ENV MAMBA_USER_ID=$MAMBA_USER_ID
ENV MAMBA_USER_GID=$MAMBA_USER_GID

COPY --from=mambaorg/micromamba:latest /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY --from=mambaorg/micromamba:latest /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh
RUN /usr/local/bin/_dockerfile_initialize_user_accounts.sh && /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN mkdir /app && chown mambauser:mambauser /app

USER mambauser

WORKDIR /app

COPY environment.yml /app/
COPY requirements.txt /app/

RUN micromamba install -n base -c conda-forge awscli -y && \
    micromamba run -n base aws s3 cp s3://datagrok-data/models/admetox /app --no-sign-request --recursive --region us-east-2 && \
    micromamba remove -n base awscli -y && \
    micromamba clean -afy

RUN set -ex ; \
    micromamba env create -f environment.yml ; \
    micromamba clean -afy

COPY . /app

SHELL ["micromamba","run","-n","myenv", "/bin/bash","-c"]

EXPOSE 8000

ENTRYPOINT ["micromamba", "run", "-n", "myenv", "gunicorn", "-w", "1", "-b", "0.0.0.0:8000", "--timeout", "3000", "ADMETox:app"]
