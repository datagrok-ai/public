# aizynthfinder Celery worker on Wolfi (glibc, near-zero CVEs). The venv is built in a builder stage
# and copied into a minimal final stage without pip/setuptools/build tools.
# onnxruntime<1.17 ships an executable-stack .so (https://github.com/microsoft/onnxruntime/issues/24911);
# patchelf clears the flag so the pre-latest aizynthfinder keeps working.

FROM cgr.dev/chainguard/wolfi-base AS builder

RUN apk add --no-cache python-3.11 py3.11-pip patchelf

COPY requirements.txt /tmp/
RUN python3.11 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt && \
    /opt/venv/bin/pip install --no-cache-dir --no-deps aizynthfinder==4.3.2 && \
    patchelf --clear-execstack /opt/venv/lib/python3.11/site-packages/onnxruntime/capi/onnxruntime_pybind11_state.cpython-311-x86_64-linux-gnu.so && \
    /opt/venv/bin/pip uninstall -y pip wheel

FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache && apk add --no-cache python-3.11 libgomp libstdc++ libxrender libxext && \
    rm -rf /usr/lib/python3.11/ensurepip && \
    adduser -D datagrok

COPY --from=builder /opt/venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH

WORKDIR /app

RUN python -c "import pyarrow, pandas, tables, onnxruntime, rdkit.Chem.Draw, datagrok_celery_task, datagrok_api" && \
    mkdir -p /app/configs/default && \
    python -c "import sys; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', '/app/configs/default']; main()" && \
    python -c "from aizynthfinder.aizynthfinder import AiZynthFinder; AiZynthFinder(configfile='/app/configs/default/config.yml')"
COPY app.py /app/

RUN chown -R datagrok:datagrok /app
USER datagrok
RUN python -c "import app"
# grok_spawner expects port to be exposed
EXPOSE 8000

# Datagrok will populate this env variables based on DockerContainer id and config properties
# Use concurrency=cpu in container.json
CMD celery -A app worker \
    --loglevel=info \
    --hostname=${CELERY_HOSTNAME:-aizynthfinder@%h} \
    --concurrency=1 \
    -Q ${TASK_QUEUE_NAME:-celery}
