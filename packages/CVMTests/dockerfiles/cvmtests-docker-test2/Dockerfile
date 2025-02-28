FROM python:3.8-alpine3.18 as pip_builder

ARG VIRTUAL_ENV=/app/venv
ENV PATH "${VIRTUAL_ENV}/bin:$PATH"

COPY requirements.txt ${VIRTUAL_ENV}/requirements.txt

RUN python -m venv ${VIRTUAL_ENV} && \
    pip install --timeout 3600 --no-cache-dir --requirement ${VIRTUAL_ENV}/requirements.txt && \
    rm -rf /tmp/*


FROM python:3.8-alpine3.18

WORKDIR /app

ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=/app/venv
RUN addgroup -g 2001 grok && \
    adduser -S -h ${HOME_DIR} -G grok -u 1001 grok && \
    chown -R grok:grok ${HOME_DIR}

COPY --chown=grok:grok --from=pip_builder $VIRTUAL_ENV $VIRTUAL_ENV
ENV PATH ${VIRTUAL_ENV}/bin:$PATH
ENV PYTHONUNBUFFERED=1
COPY . .
EXPOSE 5353
CMD ["python", "app.py"]