FROM python:3.11-slim-bookworm as pip_builder

ARG VIRTUAL_ENV=/app/venv
ENV PATH "${VIRTUAL_ENV}/bin:$PATH"

COPY requirements.in ${VIRTUAL_ENV}/requirements.in

RUN python -m venv ${VIRTUAL_ENV} && \
    pip install --timeout 3600 --no-cache-dir -r ${VIRTUAL_ENV}/requirements.in && \
    rm -rf /tmp/*


FROM python:3.11-slim-bookworm

WORKDIR /app

ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=/app/venv

COPY --from=pip_builder $VIRTUAL_ENV $VIRTUAL_ENV

ENV PATH ${VIRTUAL_ENV}/bin:$PATH
ENV CHEM_HOST='0.0.0.0'
ENV CHEM_PORT=5000
ENV CHEM_NUM_CORES=4

COPY . .
ARG CHEM_USER=chem
RUN adduser --system --group $CHEM_USER
USER $CHEM_USER
EXPOSE 5000
SHELL ["/bin/bash", "-c"]

ENTRYPOINT gunicorn app:app \
             --timeout 900 \
             -b $CHEM_HOST:$CHEM_PORT \
             --access-logfile - \
             --error-logfile - \
             --workers 2 \
             --threads=4 \
             --worker-class=gthread
