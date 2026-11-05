# Wolfi (glibc, near-zero CVEs) instead of Debian: `mlflow models serve` keeps its default
# virtualenv env manager, which needs pyenv, a C toolchain and the -dev libraries at runtime
# (pyenv builds the model's Python when it is not one of the pre-built versions, pip builds sdists).
FROM cgr.dev/chainguard/wolfi-base

ENV MLFLOW_TRACKING_URI=http://192.168.1.101:5000
ENV MODEL_NAME=datagrok-linear-regression
ENV MODEL_VERSION=1

RUN apk upgrade --no-cache && apk add --no-cache \
    bash git curl wget python-3.11 py3.11-pip \
    build-base openssl-dev zlib-dev bzip2-dev readline-dev sqlite-dev ncurses-dev xz-dev \
    libffi-dev tk-dev libxml2-dev xmlsec-dev

# cryptography capped <51 to stay within mlflow's own cryptography<51 requirement
RUN python3.11 -m venv --without-pip /opt/mlflow && \
    python3.11 -m pip --python /opt/mlflow/bin/python install --no-cache-dir \
      "mlflow[extras,azure,db,gateway,genai,auth]==3.17.2" boto3 "cryptography>=50.0.0,<51" && \
    python3.11 -m pip download --no-cache-dir --no-deps --only-binary=:all: -d /opt/python-wheels pip setuptools && \
    apk del --no-cache py3.11-pip && \
    rm -rf /usr/lib/python3.11/ensurepip
ENV PATH="/opt/mlflow/bin:/root/.pyenv/bin:${PATH}"

# Every pyenv Python (pre-built here or built by mlflow at serve time) bootstraps pip for the
# model's `python -m venv` from the current wheels in WHEEL_PKG_DIR instead of its stale
# ensurepip/_bundled copies, which are removed.
ENV PYTHON_CONFIGURE_OPTS="--with-wheel-pkg-dir=/opt/python-wheels"
RUN git clone --depth 1 https://github.com/pyenv/pyenv.git /root/.pyenv && \
    rm -rf /root/.pyenv/.git && \
    for v in 3.10.15 3.12.2 3.11.9; do pyenv install $v || exit 1; done && \
    rm -rf /root/.pyenv/versions/*/lib/python*/ensurepip/_bundled/*.whl \
      /root/.pyenv/versions/*/lib/python*/test/wheeldata/*.whl /tmp/python-build.* && \
    for v in /root/.pyenv/versions/*; do \
      $v/bin/python -m venv /tmp/venv-check && /tmp/venv-check/bin/python -m pip --version && rm -rf /tmp/venv-check || exit 1; \
    done && \
    mlflow --version

COPY run.sh run.sh
RUN chmod +x run.sh

EXPOSE 5001
CMD ["./run.sh"]
