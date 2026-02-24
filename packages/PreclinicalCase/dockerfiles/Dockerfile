FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV CORE_BINARY=/app/core
WORKDIR /app

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    unzip \
    jq \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/clinicalcase-requirements.txt

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r /tmp/clinicalcase-requirements.txt && \
    rm /tmp/clinicalcase-requirements.txt

RUN useradd -m datagrok

RUN set -euo pipefail \
    && LATEST_RELEASE_URL=$(curl -s --fail --retry 3 https://api.github.com/repos/cdisc-org/cdisc-rules-engine/releases/latest | jq -r '.assets[] | select(.name == "core-ubuntu-latest.zip") | .browser_download_url') \
    && echo "Downloading from: ${LATEST_RELEASE_URL}" \
    && curl -L \
        --fail \
        --retry 10 \
        --retry-delay 5 \
        --retry-max-time 1800 \
        --retry-connrefused \
        --retry-all-errors \
        --connect-timeout 90 \
        -C - \
        -o /tmp/core-ubuntu-latest.zip \
        "${LATEST_RELEASE_URL}" \
    && unzip /tmp/core-ubuntu-latest.zip -d /tmp/cdisc-rules-engine \
    && rm /tmp/core-ubuntu-latest.zip \
    && mv /tmp/cdisc-rules-engine/core/* /app \
    && rm -rf /tmp/cdisc-rules-engine \
    && chmod +x /app/core

COPY app.py /app/

RUN mkdir -p /data /var/log/clinicalcase && \
    chown -R datagrok:datagrok /app /data /var/log/clinicalcase

USER datagrok
EXPOSE 8000

CMD celery -A app worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=1 \
    -Q $TASK_QUEUE_NAME