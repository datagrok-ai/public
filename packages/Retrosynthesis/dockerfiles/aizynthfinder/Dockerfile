FROM python:3.10 as builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libpq5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --upgrade pip && \
    pip install --prefix=/install -r requirements.txt && \
    pip install --prefix=/install aizynthfinder[all]

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
EXPOSE 38001

# Datagrok will populate this env variables based on DockerContainer id and config properties
# Use concurrency=cpu in container.json
CMD celery -A app worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=1 \
    -Q $TASK_QUEUE_NAME
