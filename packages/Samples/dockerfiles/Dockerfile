# The conda env is built with micromamba in a builder stage and copied onto Wolfi (glibc, near-zero CVEs);
# the final image carries only the env, without micromamba, pip or the Debian base.
FROM mambaorg/micromamba:1.5.3 AS builder

COPY environment.yaml /tmp/environment.yaml

RUN micromamba create -y -p /opt/conda/envs/myenv -f /tmp/environment.yaml && \
    /opt/conda/envs/myenv/bin/python -m pip uninstall -y pip wheel && \
    rm -rf /opt/conda/envs/myenv/lib/python*/ensurepip && \
    micromamba clean --all --yes

FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache

COPY --from=builder /opt/conda/envs/myenv /opt/conda/envs/myenv
ENV PATH=/opt/conda/envs/myenv/bin:$PATH

WORKDIR /app

COPY . .

RUN python -c "import app"

EXPOSE 5555

ENTRYPOINT celery -A app worker \
    --loglevel=info \
    --hostname=${CELERY_HOSTNAME:-samples@%h} \
    --concurrency=1 \
    -Q ${TASK_QUEUE_NAME:-celery}
