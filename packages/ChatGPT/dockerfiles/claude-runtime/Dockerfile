FROM python:3.11-slim AS pip_builder

ARG VIRTUAL_ENV=/app/venv
ENV PATH="${VIRTUAL_ENV}/bin:$PATH"

COPY requirements.txt ${VIRTUAL_ENV}/requirements.txt

RUN python -m venv ${VIRTUAL_ENV} && \
    pip install --timeout 3600 --no-cache-dir --requirement ${VIRTUAL_ENV}/requirements.txt && \
    rm -rf /tmp/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl ca-certificates && \
    curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y --no-install-recommends nodejs && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

FROM python:3.11-slim

WORKDIR /app

ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=/app/venv
RUN groupadd -g 2001 grok && \
    useradd -m -d ${HOME_DIR} -g grok -u 1001 grok && \
    mkdir -p /workspace && \
    chown -R grok:grok ${HOME_DIR} /workspace

COPY --chown=grok:grok --from=pip_builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}
COPY --from=pip_builder /usr/bin/node /usr/bin/node
COPY --from=pip_builder /usr/lib/node_modules /usr/lib/node_modules
RUN ln -sf /usr/lib/node_modules/npm/bin/npm-cli.js /usr/bin/npm

ENV PATH="${VIRTUAL_ENV}/bin:$PATH"
ENV PYTHONUNBUFFERED=1

# Anthropic API key — override at runtime via -e or container env config
ENV ANTHROPIC_API_KEY=""

# Workspace directory for the mounted public folder
# Mount via: docker run -v /path/to/public:/workspace ...
ENV CLAUDE_WORKSPACE=/workspace

COPY --chown=grok:grok . .

USER grok

EXPOSE 5353

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD wget --spider http://127.0.0.1:5353/health || exit 1

CMD ["python", "app.py"]
