FROM mambaorg/micromamba:2.10.0 AS builder

ENV PIP_NO_CACHE_DIR=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# padelpy's bundled PaDEL-Descriptor Java tool (log4j 1.2.15, guava 17.0) comes in via
# chemprop -> astartes[molecules] -> aimsim -> padelpy and is unused at inference.
RUN set -ex ; \
    micromamba env create -f /tmp/environment.yml -y ; \
    env=$MAMBA_ROOT_PREFIX/envs/myenv ; \
    site=$(echo $env/lib/python3.*/site-packages) ; \
    rm -rf $site/pip $site/pip-* $site/setuptools $site/setuptools-* $site/_distutils_hack \
        $site/distutils-precedence.pth $site/wheel $site/wheel-* $site/pkg_resources \
        $site/padelpy/PaDEL-Descriptor $env/lib/python3.*/ensurepip/_bundled ; \
    find $env/conda-meta -name 'pip-*.json' -o -name 'setuptools-*.json' -o -name 'wheel-*.json' | xargs -r rm -f

RUN mkdir /tmp/models && $MAMBA_ROOT_PREFIX/envs/myenv/bin/python -c "\
import boto3; from botocore import UNSIGNED; from botocore.config import Config; \
s3 = boto3.client('s3', region_name='us-east-2', config=Config(signature_version=UNSIGNED)); \
keys = [o['Key'] for o in s3.list_objects_v2(Bucket='datagrok-data', Prefix='models/admetox/')['Contents']]; \
[s3.download_file('datagrok-data', k, '/tmp/models/' + k.rsplit('/', 1)[1]) for k in keys if not k.endswith('/')]"

FROM cgr.dev/chainguard/wolfi-base

ARG MAMBA_USER=mambauser
ARG MAMBA_USER_ID=57439
ARG MAMBA_USER_GID=57439
ENV MAMBA_USER=$MAMBA_USER \
    MAMBA_USER_ID=$MAMBA_USER_ID \
    MAMBA_USER_GID=$MAMBA_USER_GID \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/myenv/bin:$PATH \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    LD_LIBRARY_PATH=/usr/local/nvidia/lib:/usr/local/nvidia/lib64

RUN apk upgrade --no-cache && apk add --no-cache libstdc++ libgcc && \
    addgroup -g $MAMBA_USER_GID $MAMBA_USER && \
    adduser -D -u $MAMBA_USER_ID -G $MAMBA_USER $MAMBA_USER

COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /opt/conda/envs/myenv /opt/conda/envs/myenv
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /tmp/models /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

USER $MAMBA_USER
WORKDIR /app

# Load a checkpoint and predict on CPU, as the worker does, so a missing library fails the build.
RUN python -c "import padelpy, chemprop, pandas, Admetica; \
r = Admetica.predict(pandas.Series(['CCO', 'c1ccccc1O']), 'Caco2,PPBR'); print(r); \
assert r.notna().all().all()" && rm -rf /app/lightning_logs

EXPOSE 8000

ENTRYPOINT exec celery -A Admetica worker \
    --loglevel=info \
    --hostname=$CELERY_HOSTNAME \
    --concurrency=1 \
    -Q $TASK_QUEUE_NAME
