# Add ability to easily change python version
ARG PYTHON_VERSION=3.8

# Use our base image which is optimised and secured by our flow 
FROM datagrok/python:${PYTHON_VERSION}

MAINTAINER Volodymyr Dyma <vdyma@datagrok.ai>

# Use pipefail to avoid errors after piping
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG HOME_DIR=/home/grok
# It is important to run docker container as user and not as root
RUN groupadd --gid 2001 grok ; \
    useradd --system --create-home --home ${HOME_DIR} --gid grok --uid 1001 grok

# Remove unzip after unpacking
ARG DEBIAN_FRONTEND=noninteractive
RUN savedAptMark="$(apt-mark showmanual)" ; \
    apt-get update ; \
    apt-get install -y --no-install-recommends \
      unzip \
      wget \
    ; \
    apt-mark auto '.*' > /dev/null ; \
    [ -z "$savedAptMark" ] || apt-mark manual $savedAptMark ; \
	wget https://mafft.cbrc.jp/alignment/software/mafft_7.490-1_amd64.deb -O mafft.deb; \
	apt install -y ./mafft.deb; \
	rm -rf mafft.deb; \
    wget https://github.com/Merck/PepSeA/archive/refs/heads/main.zip -O PepSeA.zip; \
    unzip -q PepSeA.zip -d /opt ; \
	rm -rf PepSeA.zip; \
    chown -R grok:grok /opt/PepSeA-main ; \
    apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false ; \
    apt-get clean ; \
    rm -rf /var/lib/apt/lists/* /var/cache/* /var/tmp/* /tmp/* /opt/PepSeA-main/alignment/mafft /opt/PepSeA-main/local.env;

USER grok:grok

# Command source does not work for Docker, cause it will apply only to one layer
# The PATH works better for Docker
ARG VIRTUAL_ENV=${HOME_DIR}/venv
ENV PATH "${VIRTUAL_ENV}/bin:$PATH"
RUN mkdir -p ${VIRTUAL_ENV} ; \
    python -m venv ${VIRTUAL_ENV} ; \
    pip install --no-cache-dir --upgrade pip ; \
    pip install --no-cache-dir --upgrade wheel setuptools ; \
    pip install --timeout 3600 --no-cache-dir \
      --requirement /opt/PepSeA-main/requirements.txt ; \
    find ${VIRTUAL_ENV} -type d \( -path *tensorflow \) -prune -false -depth \( \
      \( -type d -a \( -name test -o -name tests -o -name idle_test \) \) \
          -o \( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
      \) -exec rm -rf '{}' + ; \
    rm -rf /var/cache/* /var/tmp/* /tmp/*

WORKDIR /opt/PepSeA-main

ENV ADD_HTTPS_MIDDLEWARE=""
ENV MAFFT_DIR="/usr/bin/"

# CMD can be changed easially, you need ENTRYPOINT in this case
ENTRYPOINT ["uvicorn", "--host", "0.0.0.0", "alignment.api:api"]

EXPOSE 8000
