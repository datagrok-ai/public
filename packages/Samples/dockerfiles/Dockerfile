FROM mambaorg/micromamba:1.5.3

USER root

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY environment.yaml .

RUN micromamba create -n myenv -f environment.yaml && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "myenv", "/bin/bash", "-c"]

COPY . .

EXPOSE 5555

ENTRYPOINT micromamba run -n myenv celery -A app worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=1 \
    -Q $TASK_QUEUE_NAME