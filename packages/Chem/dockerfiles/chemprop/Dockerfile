FROM cgr.dev/chainguard/wolfi-base AS pip_builder

ARG VIRTUAL_ENV=/opt/venv

RUN apk add --no-cache python-3.12 py3.12-pip libstdc++

COPY requirements.in ${VIRTUAL_ENV}/requirements.in

# CPU torch: our hosts have no GPUs, and the service and the chemprop CLI (Lightning's
# auto accelerator) fall back to CPU. Installed first so the CUDA build from PyPI is not pulled.
# padelpy's bundled PaDEL-Descriptor Java tool (log4j 1.2, guava jars) is dropped: it comes via
# chemprop -> astartes[molecules] -> aimsim -> padelpy and is unused at inference.
RUN python3.12 -m venv --without-pip ${VIRTUAL_ENV} && \
    python3.12 -m pip --python ${VIRTUAL_ENV}/bin/python install --no-cache-dir \
      --index-url https://download.pytorch.org/whl/cpu 'torch>=2.9' && \
    python3.12 -m pip --python ${VIRTUAL_ENV}/bin/python install --timeout 3600 --no-cache-dir \
      -r ${VIRTUAL_ENV}/requirements.in && \
    rm -rf ${VIRTUAL_ENV}/lib/python3.12/site-packages/padelpy/PaDEL-Descriptor && \
    find ${VIRTUAL_ENV} -depth -type d -name __pycache__ -exec rm -rf {} +


# Wolfi (glibc) with only the interpreter: the venv is copied in without pip/setuptools/wheel.
FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache && \
    apk add --no-cache python-3.12 libstdc++ libgomp && \
    rm -rf /usr/lib/python3.12/ensurepip

ARG VIRTUAL_ENV=/opt/venv
COPY --from=pip_builder $VIRTUAL_ENV $VIRTUAL_ENV

ENV PATH=${VIRTUAL_ENV}/bin:$PATH
# chemprop 2.0 checkpoints pickle classes; torch >= 2.6 defaults torch.load to weights_only=True
ENV TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1

WORKDIR /app

COPY . .
RUN python -c "import torch, lightning, padelpy, flask_cors, grok_compute" && chemprop --help > /dev/null

EXPOSE 45002

ENTRYPOINT ["python", "grok_compute.py"]
