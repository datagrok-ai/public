# node 20 is EOL since 2026-04; npm@latest (12+) needs node >=22.22
FROM node:22-trixie-slim

RUN apt-get update && apt-get -y upgrade && apt-get install -y --no-install-recommends \
    git ca-certificates python3 python3-pip util-linux && \
    python3 -m pip install --no-cache-dir --break-system-packages --ignore-installed \
      --upgrade "pip>=26.1.2" "setuptools>=78.1.1" virtualenv && \
    apt-get purge -y python3-pip python3-setuptools python3-wheel \
      python3-pip-whl python3-setuptools-whl && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    for p in pip setuptools; do \
      ls -1 /usr/local/lib/python3*/dist-packages/virtualenv/seed/wheels/embed/$p-*.whl \
        | sort -V | head -n -1 | xargs -r rm -f; \
    done

WORKDIR /app

ARG HOME_DIR=/home/grok
RUN groupadd -g 2001 grok && \
    useradd -m -d ${HOME_DIR} -g grok -u 1001 grok && \
    groupadd -g 2002 broker && \
    useradd -m -d /home/broker -g broker -u 1002 broker && \
    chmod 0700 /home/broker && \
    git clone --depth 1 https://github.com/datagrok-ai/public.git /workspace && \
    find /workspace \( -name '*.jar' -o -name 'Cargo.lock' \) -delete && \
    mkdir -p /users && \
    chown -R grok:grok ${HOME_DIR} /workspace /users

RUN dpkg --purge --force-remove-essential --force-depends \
      perl perl-base perl-modules-5.40 libperl5.40

RUN npm install -g npm@latest && \
    npm install --prefix /tmp/fix --no-save \
      ip-address@10.3.1 undici@6.28.0 brace-expansion@5.0.9 tar@7.5.21 && \
    for p in ip-address undici brace-expansion tar; do \
      rm -rf "$(npm root -g)/npm/node_modules/$p"; \
      mv "/tmp/fix/node_modules/$p" "$(npm root -g)/npm/node_modules/$p"; \
    done && \
    rm -rf /tmp/fix && \
    npm install -g @anthropic-ai/claude-code@2.1.231 && \
    node "$(npm root -g)/@anthropic-ai/claude-code/install.cjs" && \
    claude --version && \
    npm cache clean --force

COPY --chown=grok:grok package.json package-lock.json* tsconfig.json ./
COPY --chown=grok:grok src ./src
RUN npm ci && npx tsc && rm -rf src tsconfig.json && npm prune --omit=dev

COPY --chown=grok:grok plugin /app/plugin
COPY entrypoint.sh /app/entrypoint.sh
RUN chmod 0755 /app/entrypoint.sh

RUN mkdir -p /home/broker/.claude && \
    printf '{}' > /home/broker/.claude/.credentials.json && \
    chown -R broker:broker /home/broker/.claude && \
    chmod 0700 /home/broker/.claude && \
    chmod 0600 /home/broker/.claude/.credentials.json && \
    mkdir -p /home/grok/.claude && chown -R grok:grok /home/grok/.claude

ENV CLAUDE_WORKSPACE=/workspace

EXPOSE 5355

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:5355/health').then(r=>r.ok?process.exit(0):process.exit(1)).catch(()=>process.exit(1))"

ENTRYPOINT ["/app/entrypoint.sh"]
