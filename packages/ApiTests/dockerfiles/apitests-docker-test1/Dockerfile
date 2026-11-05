FROM python:3.12-alpine as pip_builder

ARG VIRTUAL_ENV=/app/venv
ENV PATH "${VIRTUAL_ENV}/bin:$PATH"

COPY requirements.txt ${VIRTUAL_ENV}/requirements.txt

RUN python -m venv ${VIRTUAL_ENV} && \
    pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --timeout 3600 --no-cache-dir --requirement ${VIRTUAL_ENV}/requirements.txt && \
    pip uninstall -y setuptools wheel pip && \
    rm -rf /tmp/*


FROM python:3.12-alpine

WORKDIR /app

# Patch base-image CVEs with apk upgrade, and drop pip/setuptools/wheel and ensurepip's
# bundled wheels: nothing installs at runtime.
RUN apk upgrade --no-cache && \
    python -m pip uninstall -y pip setuptools wheel && \
    rm -rf /usr/local/lib/python3*/ensurepip

ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=/app/venv
RUN addgroup -g 2001 grok && \
    adduser -S -h ${HOME_DIR} -G grok -u 1001 grok && \
    chown -R grok:grok ${HOME_DIR}

COPY --chown=grok:grok --from=pip_builder $VIRTUAL_ENV $VIRTUAL_ENV
ENV PATH ${VIRTUAL_ENV}/bin:$PATH

COPY . .
EXPOSE 5353
CMD ["python", "app.py"]
