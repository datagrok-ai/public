FROM python:3.10-slim-trixie as base

ENV MLFLOW_TRACKING_URI=http://192.168.1.101:5000
ENV MODEL_NAME=datagrok-linear-regression
ENV MODEL_VERSION=1

# ghcr.io/mlflow/mlflow:latest (= v3.14.0 full) is still built on EOL Debian 11 bullseye,
# so replicate its single mlflow layer on a supported Debian 13 python base instead
RUN apt-get -y update && apt-get -y upgrade
RUN pip install --no-cache-dir "mlflow[extras,azure,db,gateway,genai,auth]==3.14.0"

COPY run.sh run.sh
RUN chmod +x run.sh

RUN apt-get -y update && apt-get install -y curl git gcc make
RUN apt-get -y update && apt-get install -y --no-install-recommends make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev
RUN apt-get install -y mecab-ipadic-utf8
RUN python -v
RUN curl https://pyenv.run | bash \
	&& python -m  pip install virtualenv
# cryptography capped <49 to stay within mlflow's own cryptography<49 requirement
RUN python -m pip install --upgrade "pip>=26.1.2" "setuptools>=78.1.1" "wheel>=0.46.2" "cryptography>=48.0.1,<49"
RUN python -m pip install boto3
ENV PATH="/root/.pyenv/bin:${PATH}"
RUN pyenv install 3.10.15
RUN pyenv install 3.12.2
RUN pyenv install 3.11.9
# refresh pip/setuptools inside the pre-baked pyenv interpreters, then drop vulnerable
# wheels pip cannot upgrade: ensurepip fallbacks, cpython test fixtures, and the
# virtualenv seed wheel used only for py3.9 targets (none exist in this image)
RUN for v in /root/.pyenv/versions/*; do $v/bin/python -m pip install --upgrade "pip>=26.1.2" "setuptools>=78.1.1"; done \
	&& rm -f /root/.pyenv/versions/*/lib/python*/ensurepip/_bundled/*.whl \
		/root/.pyenv/versions/*/lib/python*/test/wheeldata/*.whl \
		/usr/local/lib/python*/ensurepip/_bundled/*.whl \
		/usr/local/lib/python*/site-packages/virtualenv/seed/wheels/embed/pip-26.0.1-py3-none-any.whl

EXPOSE 5001
CMD ["./run.sh"]
