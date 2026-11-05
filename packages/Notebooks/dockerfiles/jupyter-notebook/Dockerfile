ARG PYTHON_VERSION=3.11
ARG CONFIG_SOURCE_DIR=./src
ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=${HOME_DIR}/venv

# Wolfi (glibc, near-zero CVEs) for every stage that reaches the final image: the venv and the
# tools are built in builder stages and only their outputs are copied onto a clean wolfi-base
FROM cgr.dev/chainguard/wolfi-base AS build_base

ARG PYTHON_VERSION
RUN apk add --no-cache \
      bash build-base curl git patch xz \
      zlib-dev gd-dev protobuf-dev protoc \
      python-${PYTHON_VERSION} python-${PYTHON_VERSION}-dev py${PYTHON_VERSION}-pip

SHELL ["/bin/bash", "-e", "-o", "pipefail", "-c"]

# antiword (textract's .doc backend, with Debian's security patches), gnuplot and pandoc are not in Wolfi
FROM build_base AS tools_builder

ARG ANTIWORD_DEBIAN=0.37-18
ARG GNUPLOT_VERSION=6.0.3
ARG PANDOC_VERSION=3.13
RUN curl -fsSL https://deb.debian.org/debian/pool/main/a/antiword/antiword_0.37.orig.tar.gz -o /tmp/antiword.tgz ; \
    curl -fsSL https://deb.debian.org/debian/pool/main/a/antiword/antiword_${ANTIWORD_DEBIAN}.debian.tar.xz -o /tmp/antiword-debian.txz ; \
    curl -fsSL https://downloads.sourceforge.net/project/gnuplot/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz -o /tmp/gnuplot.tgz ; \
    curl -fsSL https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-linux-amd64.tar.gz -o /tmp/pandoc.tgz ; \
    printf '%s\n' \
      '4415e79d9f4c8d282a1cffbdaffe7ec0178982b9608e79bfd18561234a43e0cc  /tmp/antiword.tgz' \
      '1f4ea3e79742584114d9d03d7457d8945b790f4bfec5abc09b4b02c13861c933  /tmp/antiword-debian.txz' \
      'ec52e3af8c4083d4538152b3f13db47f6d29929a3f6ecec5365c834e77f251ab  /tmp/gnuplot.tgz' \
      '38aa33b325a2bc96a4c39409631212db55702f564250bcdbead180b3bf646389  /tmp/pandoc.tgz' \
      | sha256sum -c - ; \
    tar -xzf /tmp/antiword.tgz -C /tmp ; \
    tar -xJf /tmp/antiword-debian.txz -C /tmp/antiword-0.37 ; \
    cd /tmp/antiword-0.37 ; \
    for p in $(cat debian/patches/series) ; do patch -p1 < debian/patches/$p ; done ; \
    # K&R-style prototypes do not compile as C23, gcc's default since 15 \
    make -f Makefile.Linux OPT='-O2 -std=gnu17 -fpermissive' DESTDIR=/out global_install ; \
    tar -xzf /tmp/gnuplot.tgz -C /tmp ; \
    cd /tmp/gnuplot-${GNUPLOT_VERSION} ; \
    ./configure --prefix=/usr/local --with-gd --with-readline=builtin \
      --without-qt --without-x --without-cairo --without-lua --without-libcerf --disable-wxwidgets ; \
    make -j"$(nproc)" ; \
    make DESTDIR=/out install ; \
    mkdir -p /out/usr/local/bin ; \
    tar -xzf /tmp/pandoc.tgz -C /tmp ; \
    cp /tmp/pandoc-${PANDOC_VERSION}/bin/pandoc /out/usr/local/bin/ ; \
    rm -rf /out/usr/local/share/man /out/usr/local/share/info /tmp/*

# TeX Live for nbconvert's PDF export (xelatex); the tree is self-contained, so any builder will do
FROM debian:bookworm-slim AS texlive

SHELL ["/bin/bash", "-e", "-o", "pipefail", "-c"]
RUN apt-get update ; \
    apt-get install -y --no-install-recommends ca-certificates perl wget xz-utils ; \
    mkdir /tmp/itl ; \
    wget -qO- https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz | tar -xz -C /tmp/itl --strip-components=1 ; \
    printf '%s\n' \
      'selected_scheme scheme-basic' \
      'TEXDIR /usr/local/texlive' \
      'TEXMFLOCAL /usr/local/texlive/texmf-local' \
      'TEXMFSYSVAR /usr/local/texlive/texmf-var' \
      'TEXMFSYSCONFIG /usr/local/texlive/texmf-config' \
      'TEXMFHOME ~/texmf' \
      'TEXMFVAR ~/.texlive/texmf-var' \
      'TEXMFCONFIG ~/.texlive/texmf-config' \
      'binary_x86_64-linux 1' \
      'instopt_adjustpath 0' \
      'instopt_letter 0' \
      'tlpdbopt_autobackup 0' \
      'tlpdbopt_install_docfiles 0' \
      'tlpdbopt_install_srcfiles 0' \
      > /tmp/texlive.profile ; \
    /tmp/itl/install-tl --profile /tmp/texlive.profile --repository https://mirror.ctan.org/systems/texlive/tlnet ; \
    # what Debian's texlive-xetex, texlive-fonts-recommended and texlive-plain-generic bring, \
    # plus latex-extra for the packages nbconvert's latex template loads (tcolorbox, adjustbox, ...) \
    /usr/local/texlive/bin/x86_64-linux/tlmgr install \
      collection-xetex collection-latexrecommended collection-latexextra \
      collection-fontsrecommended collection-plaingeneric ; \
    rm -rf /usr/local/texlive/tlpkg/backups /usr/local/texlive/install-tl* /usr/local/texlive/texmf-dist/doc /tmp/*

FROM build_base AS requirements_builder

ARG CONFIG_SOURCE_DIR
ARG PYTHON_VERSION
COPY ${CONFIG_SOURCE_DIR}/requirements.in /tmp/requirements.in
# build-base + python headers: some pinned packages (pyphysio, image, ...)
# ship only sdists, so pip-compile must build them for metadata
RUN python${PYTHON_VERSION} -m venv --without-pip /tmp/tools ; \
    python${PYTHON_VERSION} -m pip --python /tmp/tools/bin/python install --no-cache-dir pip pip-tools ; \
    # torch_stable.html is frozen at torch 2.3.x; +cpu wheels live on the per-package cpu indexes \
    /tmp/tools/bin/pip-compile \
      --pip-args='--retries 5 --timeout 300' \
      --find-links='https://download.pytorch.org/whl/cpu/torch/' \
      --find-links='https://download.pytorch.org/whl/cpu/torchvision/' \
      /tmp/requirements.in > /tmp/requirements.txt

FROM build_base AS pip_builder

ARG VIRTUAL_ENV
ARG PYTHON_VERSION

ENV PATH="/usr/local/kyotocabinet/bin/:${VIRTUAL_ENV}/bin:$PATH"

COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/include/kc*.h /usr/local/include/
COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/lib/libkyotocabinet.* /usr/lib/
COPY --from=datagrok/kyotocabinet:1.2.79 /usr/local/bin/kc* /usr/local/kyotocabinet/bin/
COPY --chown=root:root --from=requirements_builder /tmp/requirements.txt ${VIRTUAL_ENV}/requirements.txt

# Workaround for textract: it uses the old version of six package, which is okay for us,
# although we can not add textract to requirments file
RUN python${PYTHON_VERSION} -m venv --without-pip ${VIRTUAL_ENV} ; \
    # pip stays in the venv: notebook users install packages with %pip \
    python${PYTHON_VERSION} -m pip --python ${VIRTUAL_ENV}/bin/python install --no-cache-dir --upgrade pip setuptools wheel ; \
    pip install --no-cache-dir textract==1.6.3 ; \
    pip install --timeout 3600 --no-cache-dir \
      --requirement ${VIRTUAL_ENV}/requirements.txt ; \
    # cld3 bindings ship no cp311 wheels; their sdists need protoc, pybind11 and \
    # (pycld3) re-cythonizing, since its pregenerated cython output predates py3.11; \
    # pycld3 asks for C++11, which current protobuf headers no longer compile with \
    pip install --no-cache-dir pybind11 'cython<3' ; \
    curl -fsSL https://files.pythonhosted.org/packages/source/p/pycld3/pycld3-0.22.tar.gz -o /tmp/pycld3.tgz ; \
    echo '2de0f6895342a116bce6032e6b5bea747de1dcad5d511c6583ae4343a0708dd5  /tmp/pycld3.tgz' | sha256sum -c - ; \
    tar -xzf /tmp/pycld3.tgz -C /tmp ; \
    sed -i 's/-std=c++11/-std=c++17/' /tmp/pycld3-0.22/setup.py ; \
    USE_CYTHON=1 pip install --no-cache-dir --no-build-isolation /tmp/pycld3-0.22 ; \
    pip install --no-cache-dir --no-build-isolation gcld3==3.0.13 ; \
    pip uninstall -y pybind11 cython ; \
    # kyotocabinet-python's setup.py uses `os` relying on the py<=3.8 stdlib \
    # distutils leaking it via `import *`; patch the import and build against \
    # the copied-in kyotocabinet headers (kcutilmgr is on PATH here) \
    curl -fsSL https://dbmx.net/kyotocabinet/pythonpkg/kyotocabinet-python-1.23.tar.gz -o /tmp/kc.tgz ; \
    tar -xzf /tmp/kc.tgz -C /tmp ; \
    sed -i '1i import os' /tmp/kyotocabinet-python-1.23/setup.py ; \
    pip install --no-cache-dir --no-build-isolation /tmp/kyotocabinet-python-1.23 ; \
    # nbclassic (latest 1.3.3) still bundles underscore 1.13.7 (CVE-2026-27601) and moment 2.29.4 \
    # (CVE-2026-17495); swap in the patch releases \
    curl -fsSL https://registry.npmjs.org/underscore/-/underscore-1.13.8.tgz -o /tmp/underscore.tgz ; \
    curl -fsSL https://registry.npmjs.org/moment/-/moment-2.31.0.tgz -o /tmp/moment.tgz ; \
    mkdir /tmp/underscore /tmp/moment ; \
    tar -xzf /tmp/underscore.tgz -C /tmp/underscore --strip-components=1 ; \
    tar -xzf /tmp/moment.tgz -C /tmp/moment --strip-components=1 ; \
    for d in $(find ${VIRTUAL_ENV} -type d -path '*/static/components/underscore') ; do \
      cp /tmp/underscore/package.json /tmp/underscore/underscore-min.js "$d"/ ; \
    done ; \
    for d in $(find ${VIRTUAL_ENV} -type d -path '*/static/components/moment') ; do \
      cp /tmp/moment/package.json /tmp/moment/moment.js "$d"/ ; \
      cp /tmp/moment/min/moment.min.js /tmp/moment/min/moment-with-locales.js /tmp/moment/min/moment-with-locales.min.js \
         /tmp/moment/min/locales.js /tmp/moment/min/locales.min.js "$d"/min/ ; \
    done ; \
    grep -q '"version": "1.13.8"' ${VIRTUAL_ENV}/lib/python${PYTHON_VERSION}/site-packages/nbclassic/static/components/underscore/package.json ; \
    grep -q '"version": "2.31.0"' ${VIRTUAL_ENV}/lib/python${PYTHON_VERSION}/site-packages/nbclassic/static/components/moment/package.json ; \
    # explicit model-wheel URLs: `spacy download` resolves an empty version here \
    # (compatibility lookup fails in the build env) and 404s \
    pip install --no-cache-dir \
      https://github.com/explosion/spacy-models/releases/download/en_core_web_sm-3.6.0/en_core_web_sm-3.6.0-py3-none-any.whl \
      https://github.com/explosion/spacy-models/releases/download/de_core_news_sm-3.6.0/de_core_news_sm-3.6.0-py3-none-any.whl \
      https://github.com/explosion/spacy-models/releases/download/xx_ent_wiki_sm-3.6.0/xx_ent_wiki_sm-3.6.0-py3-none-any.whl ; \
    python -c "import cld3, gcld3, kyotocabinet" ; \
    find ${VIRTUAL_ENV} -type d \( -path *tensorflow \) -prune -false -depth \( \
      \( -type d -a \( -name test -o -name tests -o -name idle_test \) \) \
          -o \( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
      \) -exec rm -rf '{}' + ; \
    rm -rf /root/.cache /tmp/*

FROM cgr.dev/chainguard/wolfi-base AS main

ENV CVM_ENVIRONMENT='docker'

ARG HOME_DIR
ARG PYTHON_VERSION

# The runtime counterparts of the Debian packages this image had; the -dev ones only ever
# served compilers, which the image does not have
RUN set -e ; \
    apk upgrade --no-cache ; \
    apk add --no-cache \
      bash procps gettext-envsubst \
      python-${PYTHON_VERSION} libstdc++ libgomp zlib libprotobuf \
      nginx-stable nginx-stable-mime-types \
      # NLP \
      tesseract tesseract-eng tesseract-osd unrtf poppler-utils \
      # NLP END \
      gdal gdal-tools ghostscript libgd \
      cairo pango poppler-glib proj gmp glib \
      libx11 libxext libxrender libsm libice libglvnd mesa-gl \
      fontconfig freetype ttf-dejavu \
    ; \
    rm -rf /usr/lib/python3*/ensurepip /etc/nginx/conf.d/* /etc/nginx/http.d ; \
    addgroup -g 2001 grok ; \
    adduser -D -h ${HOME_DIR} -G grok -u 1001 grok ; \
    addgroup -g 1000 node ; \
    adduser -D -h /home/node -G node -u 1000 -s /bin/bash node ; \
    mkdir -p ${HOME_DIR}/.local/share/jupyter/kernels ; \
    mkdir -p ${HOME_DIR}/grok_helper/cache ${HOME_DIR}/grok_helper/environments ; \
    mkdir -p ${HOME_DIR}/notebooks ; \
    mkdir -p ${HOME_DIR}/.jupyter ; \
    chown -R grok:grok ${HOME_DIR}

SHELL ["/bin/bash", "-e", "-o", "pipefail", "-c"]

ARG VIRTUAL_ENV
COPY --from=tools_builder /out/usr/local/ /usr/local/
COPY --from=texlive /usr/local/texlive /usr/local/texlive
COPY --chown=grok:grok --from=pip_builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}
COPY --chown=root:root --from=pip_builder /usr/local/include/kc*.h /usr/local/include/
COPY --chown=root:root --from=pip_builder /usr/lib/libkyotocabinet.* /usr/lib/
COPY --chown=root:root --from=pip_builder /usr/local/kyotocabinet/bin/ /usr/local/kyotocabinet/bin/

ENV PATH="${VIRTUAL_ENV}/bin:/usr/local/kyotocabinet/bin:/usr/local/texlive/bin/x86_64-linux:$PATH"

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
RUN mkdir -p /etc/nginx/sites-enabled /etc/nginx/conf.d /var/lib/nginx/tmp /var/log/nginx /var/cache/nginx /run/nginx

COPY nginx.conf.template /etc/nginx/nginx.conf.template
COPY grok_cvm.conf.template /etc/nginx/sites-enabled/grok_cvm.conf
COPY entrypoint.sh /etc/nginx/entrypoint.sh
RUN sed -i 's/\r$//' /etc/nginx/entrypoint.sh && chmod +x /etc/nginx/entrypoint.sh && \
    chown -R grok:grok /etc/nginx /var/lib/nginx /var/log/nginx /var/cache/nginx /run/nginx

EXPOSE 8090

USER grok:grok

# What the image promises, checked where it runs: the tools, nginx's config, and a notebook
# executed by the kernel and exported to PDF exactly as grok_helper's /notebooks/convert does
RUN cd /tmp ; \
    command -v antiword unrtf pdftotext xelatex ; \
    tesseract --list-langs ; gdalinfo --version ; gs --version ; pandoc --version > /dev/null ; \
    gnuplot -e 'set terminal dumb; plot sin(x)' > /dev/null ; \
    RESOLVER_ADDR=127.0.0.1 envsubst '${RESOLVER_ADDR}' < /etc/nginx/nginx.conf.template > /tmp/nginx.conf ; \
    nginx -t -c /tmp/nginx.conf ; \
    python -c "import os, nbformat; from nbconvert import PDFExporter; from nbconvert.preprocessors import ExecutePreprocessor; v4 = nbformat.v4; nb = v4.new_notebook(cells=[v4.new_markdown_cell('# Smoke *test*'), v4.new_code_cell('import numpy, pandas, matplotlib, cv2, rdkit.Chem, sklearn, kyotocabinet, cld3, gcld3, textract, flair, torch, tensorflow, spacy; spacy.load(\"en_core_web_sm\"); print(numpy.arange(3))')]); ExecutePreprocessor(timeout=600, kernel_name='python3').preprocess(nb, {'metadata': {'path': '/tmp'}}); pdf, _ = PDFExporter(extra_template_basedirs=[os.environ['GROK_HELPER_DIR'] + '/nbconvert_templates'], template_name='hide-input-pdf').from_notebook_node(nb); assert pdf.startswith(b'%PDF'), pdf[:200]" ; \
    rm -rf /tmp/* /tmp/.[!.]* ${HOME_DIR}/.cache ${HOME_DIR}/.local/share/jupyter/runtime ${HOME_DIR}/.ipython ${HOME_DIR}/.keras

WORKDIR '/etc/nginx'
ENTRYPOINT ["/bin/sh", "/etc/nginx/entrypoint.sh"]
