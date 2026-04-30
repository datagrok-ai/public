FROM node:20-slim

RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates python3 make g++ && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

ARG HOME_DIR=/home/grok
RUN groupadd -g 2001 grok && \
    useradd -m -d ${HOME_DIR} -g grok -u 1001 grok && \
    git clone --depth 1 https://github.com/datagrok-ai/public.git /workspace && \
    mkdir -p /users && \
    chown -R grok:grok ${HOME_DIR} /workspace /users

COPY --chown=grok:grok package.json package-lock.json* tsconfig.json ./
COPY --chown=grok:grok src ./src
RUN npm ci && npx tsc && rm -rf src tsconfig.json && npm prune --production && \
    npm install -g @zilliz/claude-context-mcp@latest

COPY --chown=grok:grok plugin /app/plugin

ENV CLAUDE_WORKSPACE=/workspace

USER grok

EXPOSE 5355

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:5355/health').then(r=>r.ok?process.exit(0):process.exit(1)).catch(()=>process.exit(1))"

CMD ["node", "dist/server.js"]
