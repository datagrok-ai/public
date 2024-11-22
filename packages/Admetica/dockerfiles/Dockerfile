FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

SHELL ["/bin/bash","-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && \
    apt install -y --no-install-recommends \
        curl bzip2 ca-certificates ; \
     apt clean

ENV MAMBA_EXE=/bin/micromamba
ARG MAMBA_VERSION=1.5.3
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
ENV MAMBA_USER_GID=$MAMBA_USER_GID
ENV MAMBA_USER_ID=$MAMBA_USER_ID
ENV MAMBA_ROOT_PREFIX=/opt/conda
COPY --from=mambaorg/micromamba:latest /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY --from=mambaorg/micromamba:latest /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh
RUN /usr/local/bin/_dockerfile_initialize_user_accounts.sh && /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN mkdir /app && chown $MAMBA_USER:$MAMBA_USER /app

USER $MAMBA_USER

WORKDIR /app

COPY environment.yml /app/

RUN micromamba install -n base -c conda-forge awscli -y && \
    micromamba run -n base aws s3 cp s3://datagrok-data/models/admetox /app --no-sign-request --recursive --region us-east-2 && \
    micromamba remove -n base awscli -y && \
    micromamba clean -afy

RUN set -ex ; \
    micromamba env create -f environment.yml -y ; \
    micromamba clean -afy

COPY . /app

SHELL ["micromamba","run","-n","myenv", "/bin/bash","-c"]

EXPOSE 8000

ENTRYPOINT ["micromamba", "run", "-n", "myenv", "gunicorn", "-w", "1", "-b", "0.0.0.0:8000", "--timeout", "3000", "Admetica:app"]