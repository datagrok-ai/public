FROM node:20-slim

WORKDIR /app

ARG HOME_DIR=/home/grok
RUN groupadd -g 2001 grok && \
    useradd -m -d ${HOME_DIR} -g grok -u 1001 grok && \
    mkdir -p /workspace && \
    chown -R grok:grok ${HOME_DIR} /workspace

COPY --chown=grok:grok package.json package-lock.json* tsconfig.json ./
COPY --chown=grok:grok src ./src
RUN npm ci && npx tsc && rm -rf src tsconfig.json && npm prune --production

ENV CLAUDE_WORKSPACE=/workspace

USER grok

EXPOSE 5353

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:5353/health').then(r=>r.ok?process.exit(0):process.exit(1)).catch(()=>process.exit(1))"

CMD ["node", "dist/server.js"]
