# Version 1:
# - We are using the pre-latest version of aizynthfinder because the latest release
#   fails to build. The issue arises from CGRTools being incompatible with PEP 517,
#   resulting in build errors.
#   Issue reference: https://github.com/cimm-kzn/CGRtools/issues/11
#
# - We are not using the Python slim image because the onnxruntime package fails with executable stack error

FROM mambaorg/micromamba:1.5.7

ARG DEBIAN_FRONTEND=noninteractive
ENV CONDA_ENV=myenv
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/$CONDA_ENV/bin:$PATH

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 libxrender1 libxext6 libsm6 libgl1 \
    build-essential python3-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN micromamba create -y -n $CONDA_ENV -c conda-forge python=3.10

COPY requirements.txt /app/

RUN micromamba run -n $CONDA_ENV pip install --upgrade pip && \
    micromamba run -n $CONDA_ENV pip install -r /app/requirements.txt && \
    micromamba run -n $CONDA_ENV pip install "aizynthfinder[all]==4.3.2"

WORKDIR /app

RUN mkdir -p /app/configs/default && \
    micromamba run -n $CONDA_ENV python -c "import sys; import aizynthfinder; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', '/app/configs/default']; main()"

COPY app.py /app/

RUN useradd -m datagrok && chown -R datagrok:datagrok /app
USER datagrok

EXPOSE 8000

ENTRYPOINT micromamba run -n myenv celery -A app worker --loglevel=info --hostname=$CELERY_HOSTNAME --concurrency=1 -Q $TASK_QUEUE_NAME

# Version 2:
# - We are using aizynthfinder version 4.0.0, released on December 5th, 2023 
#   (almost 2 years ago). This version does not depend on CGRTools and uses an 
#   older version of onnxruntime, which can be built successfully on slim Python images.

FROM python:3.10 as builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --upgrade pip && \
    pip install --prefix=/install -r requirements.txt && \
    pip install --prefix=/install "aizynthfinder[all]==4.0.0"

FROM python:3.10-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN useradd -m datagrok

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 libxrender1 libxext6 libsm6 libgl1 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /install/lib /usr/local/lib
COPY --from=builder /install/bin /usr/local/bin
RUN mkdir -p /app/configs/default && \
    python -c "import sys; import aizynthfinder; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', '/app/configs/default']; main()"
COPY app.py /app/

RUN chown -R datagrok:datagrok /app
USER datagrok
# grok_spawner expects port to be exposed
EXPOSE 8000

# Datagrok will populate this env variables based on DockerContainer id and config properties
# Use concurrency=cpu in container.json
CMD celery -A app worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=1 \
    -Q $TASK_QUEUE_NAME

# Version 3
# Apply the workaround from the Microsoft ONNX Runtime issue 
# (https://github.com/microsoft/onnxruntime/issues/24911) to fix the executable stack error.
# This lets us continue using a pre-latest version of aizynthfinder while keeping the image smaller 
# than a micromamba-based build. Note: this is just a workaround.

FROM python:3.10 as builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --upgrade pip && \
    pip install --prefix=/install -r requirements.txt && \
    pip install --prefix=/install "aizynthfinder[all]==4.3.2"

FROM python:3.10-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN useradd -m datagrok

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 libxrender1 libxext6 libsm6 libgl1 patchelf \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /install/lib /usr/local/lib
COPY --from=builder /install/bin /usr/local/bin
RUN patchelf --clear-execstack /usr/local/lib/python3.10/site-packages/onnxruntime/capi/onnxruntime_pybind11_state.cpython-310-x86_64-linux-gnu.so
RUN mkdir -p /app/configs/default && \
    python -c "import sys; import aizynthfinder; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', '/app/configs/default']; main()"
COPY app.py /app/

RUN chown -R datagrok:datagrok /app
USER datagrok
# grok_spawner expects port to be exposed
EXPOSE 8000

# Datagrok will populate this env variables based on DockerContainer id and config properties
# Use concurrency=cpu in container.json
CMD celery -A app worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=1 \
    -Q $TASK_QUEUE_NAME