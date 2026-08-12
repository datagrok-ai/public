ARG PYTHON_VERSION=3.11
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
# build-essential + python headers: some pinned packages (pyphysio, image, ...)
# ship only sdists, so pip-compile must build them for metadata
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
    # torch_stable.html is frozen at torch 2.3.x; +cpu wheels live on the per-package cpu indexes \
    pip-compile \
      --pip-args='--retries 5 --timeout 300' \
      --find-links='https://download.pytorch.org/whl/cpu/torch/' \
      --find-links='https://download.pytorch.org/whl/cpu/torchvision/' \
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
      protobuf-compiler \
      libprotobuf-dev \
      python${PYTHON_VERSION}-dev \
      libpython${PYTHON_VERSION}-dev \
    ; \
    python -m venv ${VIRTUAL_ENV} ; \
    pip install --no-cache-dir --upgrade pip ; \
    pip install --no-cache-dir --upgrade wheel setuptools ; \
    pip install --no-cache-dir textract==1.6.3 ; \
    pip install --timeout 3600 --no-cache-dir \
      --requirement ${VIRTUAL_ENV}/requirements.txt ; \
    # cld3 bindings ship no cp311 wheels; their sdists need protoc, pybind11 and \
    # (pycld3) re-cythonizing, since its pregenerated cython output predates py3.11 \
    pip install --no-cache-dir pybind11 'cython<3' ; \
    USE_CYTHON=1 pip install --no-cache-dir --no-build-isolation pycld3==0.22 ; \
    pip install --no-cache-dir --no-build-isolation gcld3==3.0.13 ; \
    pip uninstall -y pybind11 cython ; \
    # kyotocabinet-python's setup.py uses `os` relying on the py<=3.8 stdlib \
    # distutils leaking it via `import *`; patch the import and build against \
    # the copied-in kyotocabinet headers (kcutilmgr is on PATH here) \
    python -c "import urllib.request; urllib.request.urlretrieve('https://dbmx.net/kyotocabinet/pythonpkg/kyotocabinet-python-1.23.tar.gz', '/tmp/kc.tgz')" ; \
    tar -xzf /tmp/kc.tgz -C /tmp ; \
    sed -i '1i import os' /tmp/kyotocabinet-python-1.23/setup.py ; \
    pip install --no-cache-dir --no-build-isolation /tmp/kyotocabinet-python-1.23 ; \
    rm -rf /tmp/kc.tgz /tmp/kyotocabinet-python-1.23 ; \
    # explicit model-wheel URLs: `spacy download` resolves an empty version here \
    # (compatibility lookup fails in the build env) and 404s \
    pip install --no-cache-dir \
      https://github.com/explosion/spacy-models/releases/download/en_core_web_sm-3.6.0/en_core_web_sm-3.6.0-py3-none-any.whl \
      https://github.com/explosion/spacy-models/releases/download/de_core_news_sm-3.6.0/de_core_news_sm-3.6.0-py3-none-any.whl \
      https://github.com/explosion/spacy-models/releases/download/xx_ent_wiki_sm-3.6.0/xx_ent_wiki_sm-3.6.0-py3-none-any.whl ; \
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
    apt-get -y upgrade ; \
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
      # gnuplot-nox: the qt-enabled default drags in vulnerable qtbase libs \
      gnuplot-nox \
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
    # security floors for the base interpreter's tooling (VEX); ensurepip's bundled \
    # pip/setuptools wheels are scanner-flagged and only seed venvs, which are \
    # created in the build stages only \
    python -m pip install --no-cache-dir --upgrade pip setuptools wheel ; \
    rm -f /usr/lib/python3*/ensurepip/_bundled/*.whl ; \
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

EXPOSE 8090

USER grok:grok
WORKDIR '/etc/nginx'
ENTRYPOINT ["/bin/sh", "/etc/nginx/entrypoint.sh"]
