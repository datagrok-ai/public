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

COPY . .
EXPOSE 5000
CMD ["python", "grok_compute.py"]