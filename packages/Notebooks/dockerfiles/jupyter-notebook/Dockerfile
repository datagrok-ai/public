ARG PYTHON_VERSION=3.8
ARG CONFIG_SOURCE_DIR=./src
ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=${HOME_DIR}/venv

FROM datagrok/python:${PYTHON_VERSION} as requirements_builder

COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/include/kc*.h /usr/local/include/
COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/lib/libkyotocabinet.* /usr/local/lib/
COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/bin/kc* /usr/local/kyotocabinet/bin/

ARG CONFIG_SOURCE_DIR
ARG PYTHON_VERSION
COPY ${CONFIG_SOURCE_DIR}/requirements.in /tmp/requirements.in
# build-essential + python headers: newer releases of pinned packages' transitive deps
# (cymem, murmurhash, ...) no longer ship cp38 wheels, so pip-compile must build sdists
RUN set -ex ; \
    savedAptMark="$(apt-mark showmanual)" ; \
    apt-get update ; \
    apt-get install -y --no-install-recommends \
      git \
      zlib1g-dev \
      build-essential \
      python${PYTHON_VERSION}-dev \
      libpython${PYTHON_VERSION}-dev \
    ; \
    pip install --no-cache-dir pip-tools ; \
    pip-compile \
      --pip-args='--retries 5 --timeout 300' \
      --find-links='https://download.pytorch.org/whl/torch_stable.html' \
      /tmp/requirements.in > /tmp/requirements.txt ; \
    pip uninstall -y pip-tools ; \
    apt-mark auto '.*' > /dev/null ; \
    [ -z "$savedAptMark" ] || apt-mark manual $savedAptMark ; \
    apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false ; \
    apt-get clean ; \
    rm -rf /var/lib/apt/lists/* /var/cache/* /var/tmp/*

FROM datagrok/python:${PYTHON_VERSION} as pip_builder

ARG VIRTUAL_ENV

ENV PATH "/usr/local/kyotocabinet/bin/:${VIRTUAL_ENV}/bin:$PATH"

COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/include/kc*.h /usr/local/include/
COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/lib/libkyotocabinet.* /usr/local/lib/
COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/bin/kc* /usr/local/kyotocabinet/bin/
COPY --chown=root:root --from=requirements_builder /tmp/requirements.txt ${VIRTUAL_ENV}/requirements.txt

# Workaround for textract: it uses the old version of six package, which is okay for us,
# although we can not add textract to requirments file
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHON_VERSION
RUN set -ex ; \
    savedAptMark="$(apt-mark showmanual)" ; \
    apt-get update ; \
    apt-get install -y --no-install-recommends \
      git \
      zlib1g-dev \
      build-essential \
      python${PYTHON_VERSION}-dev \
      libpython${PYTHON_VERSION}-dev \
    ; \
    python -m venv ${VIRTUAL_ENV} ; \
    pip install --no-cache-dir --upgrade pip ; \
    pip install --no-cache-dir --upgrade wheel setuptools ; \
    pip install --no-cache-dir textract==1.6.3 ; \
    pip install --timeout 3600 --no-cache-dir \
      --requirement ${VIRTUAL_ENV}/requirements.txt ; \
    python -m spacy download en_core_web_sm ; \
    python -m spacy download de_core_news_sm ; \
    python -m spacy download xx_ent_wiki_sm ; \
    find ${VIRTUAL_ENV} -type d \( -path *tensorflow \) -prune -false -depth \( \
      \( -type d -a \( -name test -o -name tests -o -name idle_test \) \) \
          -o \( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
      \) -exec rm -rf '{}' + ; \
    apt-mark auto '.*' > /dev/null ; \
    [ -z "$savedAptMark" ] || apt-mark manual $savedAptMark ; \
    apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false ; \
    apt-get clean ; \
    rm -rf /var/lib/apt/lists/* /var/cache/* /var/tmp/* /tmp/*

FROM datagrok/python:${PYTHON_VERSION} as main

ENV CVM_ENVIRONMENT 'docker'

# We set the default shell to bash with pipefail flag called
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG HOME_DIR

ENV DEBIAN_FRONTEND=noninteractive
RUN set -ex ; \
    wget https://github.com/jgm/pandoc/releases/download/3.1.6/pandoc-3.1.6-1-amd64.deb ; \
    dpkg -i pandoc-3.1.6-1-amd64.deb ; \
    apt-get update ; \
    apt-get install -y --no-install-recommends \
      texlive-xetex texlive-fonts-recommended texlive-plain-generic \
      libprotobuf-dev \
      protobuf-compiler \
      # NLP \
      antiword \
      libtesseract4 \
      tesseract-ocr \
      unrtf \
      libx11-dev \
      # NLP END \
      gdal-bin \
      ghostscript \
      gnuplot \
      libcairo2-dev \
      libglu1-mesa-dev \
      libgdal-dev \
      libgmp-dev \
      libpango1.0-0 \
      libproj-dev \
      libpoppler-glib-dev \
      libtinfo-dev \
      libglib2.0-0 \
      zlib1g-dev \
      libgl1-mesa-glx \
    ; \
    apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false ; \
    apt-get clean ; \
    rm -rf /var/lib/apt/lists/* /var/cache/* /var/tmp/* /tmp/* ; \
    groupadd --gid 2001 grok ; \
    useradd --system --create-home --home ${HOME_DIR} --gid grok --uid 1001 grok ; \
    groupadd --gid 1000 node ; \
    useradd --uid 1000 --gid node --shell /bin/bash --create-home node ; \
    mkdir -p ${HOME_DIR}/.local/share/jupyter/kernels ; \
    mkdir -p ${HOME_DIR}/grok_helper/cache ${HOME_DIR}/grok_helper/environments ; \
    mkdir -p ${HOME_DIR}/notebooks ; \
    mkdir -p ${HOME_DIR}/.jupyter ; \
    chown -R grok:grok ${HOME_DIR}

ARG VIRTUAL_ENV
COPY --chown=grok:grok --from=pip_builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}
COPY --chown=root:root --from=pip_builder /usr/local/include/kc*.h /usr/local/include/
COPY --chown=root:root --from=pip_builder /usr/local/lib/libkyotocabinet.* /usr/local/lib/
COPY --chown=root:root --from=pip_builder /usr/local/kyotocabinet/bin/ /usr/local/kyotocabinet/bin/

ENV PATH ${VIRTUAL_ENV}/bin:/usr/local/kyotocabinet/bin:$PATH

RUN jupyter contrib nbextension install --system ; \
    jupyter nbextension enable hide_input/main ; \
    rm -rf /var/cache/* /var/tmp/* /tmp/*

ARG CONFIG_SOURCE_DIR
COPY --chown=grok:grok jupyter_notebook_config.py ${HOME_DIR}/.jupyter/
COPY --chown=grok:grok $CONFIG_SOURCE_DIR/custom.js ${HOME_DIR}/.jupyter/custom/

# GrokHelper Server
ENV GROK_HELPER_DIR=${HOME_DIR}/grok_helper
COPY --chown=grok:grok grok_helper ${GROK_HELPER_DIR}
COPY --chown=grok:grok grok_helper_configuration.yaml \
                       ${GROK_HELPER_DIR}/grok_helper_configuration.yaml
ENV GROK_HELPER_CONFIGURATION=${GROK_HELPER_DIR}/grok_helper_configuration.yaml

COPY --chown=grok:grok entrypoint-multiple-process.sh ${HOME_DIR}/entrypoint.sh
RUN sed -i 's/\r$//' ${HOME_DIR}/entrypoint.sh && chmod +x ${HOME_DIR}/entrypoint.sh

# Nginx proxy in front of jupyter + grok_helper
RUN apt-get update && apt-get install -y --no-install-recommends nginx gettext && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    mkdir -p /etc/nginx/sites-enabled /etc/nginx/conf.d

COPY nginx.conf.template /etc/nginx/nginx.conf.template
COPY grok_cvm.conf.template /etc/nginx/sites-enabled/grok_cvm.conf
COPY entrypoint.sh /etc/nginx/entrypoint.sh
RUN sed -i 's/\r$//' /etc/nginx/entrypoint.sh && chmod +x /etc/nginx/entrypoint.sh && \
    chown -R grok:grok /etc/nginx /var/lib/nginx /var/log/nginx

EXPOSE 8889 5005 8090

USER grok:grok
WORKDIR '/etc/nginx'
ENTRYPOINT ["/bin/sh", "/etc/nginx/entrypoint.sh"]
