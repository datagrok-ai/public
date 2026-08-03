FROM python:3.11-slim-bookworm as pip_builder

ARG VIRTUAL_ENV=/app/venv
ENV PATH "${VIRTUAL_ENV}/bin:$PATH"

COPY requirements.in ${VIRTUAL_ENV}/requirements.in

RUN python -m venv ${VIRTUAL_ENV} && \
    pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --timeout 3600 --no-cache-dir -r ${VIRTUAL_ENV}/requirements.in && \
    rm -rf /tmp/*


FROM python:3.11-slim-bookworm

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
# Base-image system tooling ships vulnerable pip/setuptools/wheel (incl. vendored
# jaraco.context/wheel); the app runs from the venv, but patch these so the image scans clean.
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

ARG CHEM_USER=chem
RUN adduser --system --group $CHEM_USER
USER $CHEM_USER

WORKDIR /app

ARG VIRTUAL_ENV=/app/venv

COPY --from=pip_builder --chown=$CHEM_USER:$CHEM_USER $VIRTUAL_ENV $VIRTUAL_ENV

ENV PATH ${VIRTUAL_ENV}/bin:$PATH
ENV CHEM_HOST='0.0.0.0'
ENV CHEM_PORT=45001
ENV CHEM_NUM_CORES=4

COPY --chown=$CHEM_USER:$CHEM_USER . .
EXPOSE 45001
SHELL ["/bin/bash", "-c"]

ENTRYPOINT gunicorn app:app \
             --timeout 900 \
             -b $CHEM_HOST:$CHEM_PORT \
             --access-logfile - \
             --error-logfile - \
             --workers 2 \
             --threads=4 \
             --worker-class=gthread
