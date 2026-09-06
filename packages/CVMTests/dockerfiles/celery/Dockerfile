FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libc-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# uv from PyPI, not a ghcr.io cross-registry COPY: every other image in this
# package builds from Docker Hub plus PyPI and succeeds on the CI stand, while
# this one never produced an image (builds #120-#124).
RUN pip install --no-cache-dir uv

RUN useradd -m datagrok

WORKDIR /app
COPY . /app
ENV PYTHONPATH="/app:${PYTHONPATH}"
RUN chown -R datagrok:datagrok /app

# Install Python dependencies using uv
RUN uv pip install --system celery datagrok-celery-task pyyaml && \
    if [ -f requirements.in ]; then uv pip install --system --no-cache -r requirements.in; \
    elif [ -f requirements.txt ]; then uv pip install --system --no-cache -r requirements.txt; fi

USER datagrok

EXPOSE 8000

CMD celery -A $DATAGROK_CELERY_NAME worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=$WORKERS_PROCESSES \
    -Q $TASK_QUEUE_NAME
