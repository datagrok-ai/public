# Wolfi (glibc, near-zero CVEs): the CDISC CORE binary is a glibc PyInstaller bundle.
FROM cgr.dev/chainguard/wolfi-base AS build

RUN apk add --no-cache bash python-3.12 py3.12-pip curl unzip jq
SHELL ["/bin/bash", "-c"]

COPY requirements.txt /tmp/clinicalcase-requirements.txt
RUN python3.12 -m venv --without-pip /opt/venv && \
    python3.12 -m pip --python /opt/venv/bin/python install --no-cache-dir -r /tmp/clinicalcase-requirements.txt

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
    && unzip -q /tmp/core-ubuntu-latest.zip -d /tmp/cdisc-rules-engine \
    && mkdir -p /app \
    && mv /tmp/cdisc-rules-engine/core/* /app \
    && chmod +x /app/core

FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache && apk add --no-cache python-3.12 && \
    rm -rf /usr/lib/python3.12/ensurepip && \
    adduser -D -u 1000 datagrok

ENV CORE_BINARY=/app/core
ENV PATH=/opt/venv/bin:$PATH
WORKDIR /app

COPY --from=build /opt/venv /opt/venv
COPY --from=build /app /app
COPY app.py /app/

RUN mkdir -p /data /var/log/clinicalcase && \
    chown -R datagrok:datagrok /app /data /var/log/clinicalcase

USER datagrok
RUN ls /app; /app/core version && (cd / && python -c "import celery.app") && python -c "import celery.app" && python -c "import app"

EXPOSE 8000

CMD exec celery -A app worker \
    --loglevel=info \
    --hostname=${CELERY_HOSTNAME:-clinicalcase@%h} \
    --concurrency=1 \
    -Q ${TASK_QUEUE_NAME:-clinicalcase}
