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