# node 20 is EOL since 2026-04; 22-slim is still bookworm and npm@latest (12+) needs node >=22.22
FROM node:22-slim

# apt-get upgrade first: the node-slim base ships whatever gnutls28/dpkg were
# current at ITS build time, not ours (gnutls28 3.7.9-2+deb12u7 fixes CRITICAL CVE-2026-42010).
# pip/setuptools: install fixed versions into /usr/local, then purge the apt copies
# so the vulnerable 23.0.1/66.1.1 dist-info no longer ships (python3-venv still
# seeds venvs from the *-whl packages).
RUN apt-get update && apt-get -y upgrade && apt-get install -y --no-install-recommends \
    git ca-certificates python3 python3-venv python3-pip make g++ && \
    python3 -m pip install --no-cache-dir --break-system-packages --upgrade "pip>=26.1.2" "setuptools>=78.1.1" && \
    apt-get purge -y python3-pip python3-setuptools python3-wheel && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

ARG HOME_DIR=/home/grok
# jar-delete: the checked-in JDBC driver jars (connectors/grok_connect/lib) carry
# Java CVEs and there is no JRE in this image to use them.
RUN groupadd -g 2001 grok && \
    useradd -m -d ${HOME_DIR} -g grok -u 1001 grok && \
    git clone --depth 1 https://github.com/datagrok-ai/public.git /workspace && \
    find /workspace -name '*.jar' -delete && \
    mkdir -p /users && \
    chown -R grok:grok ${HOME_DIR} /workspace /users

# Refresh npm first: the node-image npm bundles tar/glob/sigstore/cross-spawn that
# get CVE-flagged; npm 12+ ships a clean tree. npm is a runtime dependency here
# (the agent shells out to it in the workspace), so it cannot be removed like elsewhere.
RUN npm install -g npm@latest && \
    npm install -g @anthropic-ai/claude-code && \
    node "$(npm root -g)/@anthropic-ai/claude-code/install.cjs" && \
    claude --version && \
    npm cache clean --force

COPY --chown=grok:grok package.json package-lock.json* tsconfig.json ./
COPY --chown=grok:grok src ./src
RUN npm ci && npx tsc && rm -rf src tsconfig.json && npm prune --omit=dev

COPY --chown=grok:grok plugin /app/plugin

# Seed an empty subscription-credentials file so the path exists and the server's
# existence check has a target. Real subscription creds are provided at runtime (a
# read-only bind mount over this file, or the in-app `claude auth login` flow, which
# rewrites it); the API-key / Bedrock / Foundry env path does not use this file.
RUN mkdir -p /home/grok/.claude && \
    printf '{}' > /home/grok/.claude/.credentials.json && \
    chown -R grok:grok /home/grok/.claude

ENV CLAUDE_WORKSPACE=/workspace

USER grok

EXPOSE 5355

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:5355/health').then(r=>r.ok?process.exit(0):process.exit(1)).catch(()=>process.exit(1))"

CMD ["node", "dist/server.js"]
