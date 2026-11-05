# Version 3
# Apply the workaround from the Microsoft ONNX Runtime issue 
# (https://github.com/microsoft/onnxruntime/issues/24911) to fix the executable stack error.
# This lets us continue using a pre-latest version of aizynthfinder while keeping the image smaller 
# than a micromamba-based build. Note: this is just a workaround.

FROM python:3.11 as builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --upgrade pip setuptools wheel && \
    pip install --prefix=/install -r requirements.txt && \
    pip install --prefix=/install --no-deps aizynthfinder==4.3.2

FROM python:3.11-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN useradd -m datagrok

RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 libxrender1 libxext6 libsm6 libgl1 patchelf \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /install/lib /usr/local/lib
COPY --from=builder /install/bin /usr/local/bin

# pip is not used at runtime, and the setuptools copy vendored in it is flagged
RUN pip install --no-cache-dir --upgrade setuptools wheel && pip uninstall -y pip
RUN patchelf --clear-execstack /usr/local/lib/python3.11/site-packages/onnxruntime/capi/onnxruntime_pybind11_state.cpython-311-x86_64-linux-gnu.so
RUN mkdir -p /app/configs/default && \
    python -c "import sys; import aizynthfinder; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', '/app/configs/default']; main()" && \
    python -c "from aizynthfinder.aizynthfinder import AiZynthFinder; AiZynthFinder(configfile='/app/configs/default/config.yml')"
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