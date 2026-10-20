ARG HERMES_COMMIT=1f6c02c07c697026888f1c66e634534fa2aae358
ARG NODE_IMAGE=node:22-trixie-slim@sha256:154ba2f4d6fec323d28e4f4bb86bba4677f1223391a1979cf521304e03a98dfa

# Hermes' own package manager builds the venv, as Hermes' Dockerfile does: CPython and uv come from
# pm/lock.json (sha256-verified) into HERMES_RUNTIME_DIR, which the venv links to, so the path is kept.
# uv is needed for the build only. Of the unstripped CPython the lock pins, only libpython is stripped: GNU strip
# breaks the x86-64 python executable, whose .dynstr lies outside any segment.
# The sealed venv records Hermes' editable install but not its .pth; hermes_project.pth puts /opt/hermes on
# the path so the agent imports Hermes from any cwd.
FROM ${NODE_IMAGE} AS hermes
ARG HERMES_COMMIT
ENV HERMES_RUNTIME_DIR=/opt/hermes/tools
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates python3 binutils && \
    git clone --filter=blob:none https://github.com/NousResearch/hermes-agent.git /opt/hermes && \
    git -C /opt/hermes checkout "$HERMES_COMMIT" && \
    cd /opt/hermes && \
    python3 -c 'from pm import ensure; ensure("python", explicit=True)' && \
    python3 -m pm.build_env --source /opt/hermes \
      --python "$(python3 -c 'from pm import installed_package; print(installed_package("python").binary)')" \
      --out /opt/hermes/.venv --sealed --extra anthropic --extra mcp && \
    echo /opt/hermes > "$(/opt/hermes/.venv/bin/python -c 'import site; print(site.getsitepackages()[0])')/hermes_project.pth" && \
    rm -rf "$HERMES_RUNTIME_DIR"/uv-* && \
    strip --strip-debug "$HERMES_RUNTIME_DIR"/python-*/lib/libpython3.*.so.1.0 && \
    git rev-parse HEAD > COMMIT && \
    rm -rf .git tests website apps web ui-tui optional-skills && \
    /opt/hermes/.venv/bin/python -c "import run_agent, anthropic, mcp, jsonschema, ruamel.yaml"

# node 20 is EOL since 2026-04; npm@latest (12+) needs node >=22.22
FROM ${NODE_IMAGE}
ARG HERMES_COMMIT
LABEL ai.datagrok.hermes.commit=$HERMES_COMMIT

RUN apt-get update && apt-get -y upgrade && apt-get install -y --no-install-recommends \
    git ca-certificates python3 python3-pip ripgrep tini && \
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
    git clone --depth 1 --filter=blob:none --sparse https://github.com/datagrok-ai/public.git /workspace && \
    git -C /workspace sparse-checkout set --no-cone '/*' '!/docusaurus/' \
      '!*.gif' '!*.png' '!*.PNG' '!*.jpg' '!*.jpeg' '!*.bmp' '!*.mp4' '!*.jar' && \
    git clone --depth 1 --filter=blob:none --sparse https://github.com/datagrok-ai/help.git /workspace/help && \
    git -C /workspace/help sparse-checkout set --no-cone '/*' '!/uploads/' \
      '!*.gif' '!*.png' '!*.PNG' '!*.jpg' '!*.jpeg' '!*.bmp' '!*.mp4' && \
    mkdir -p /users /var/lib/grokky && \
    chown -R grok:grok ${HOME_DIR} /workspace /users /var/lib/grokky

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
    npm cache clean --force

COPY package.json package-lock.json* tsconfig.json ./
COPY src ./src
RUN npm ci && npx tsc && rm -rf src tsconfig.json && npm prune --omit=dev

COPY --from=hermes /opt/hermes /opt/hermes

COPY plugin /app/plugin
COPY agent /app/agent
COPY entrypoint.sh /app/entrypoint.sh
RUN chmod 0755 /app/entrypoint.sh

ENV CLAUDE_WORKSPACE=/workspace

EXPOSE 5355

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:5355/health').then(r=>r.ok?process.exit(0):process.exit(1)).catch(()=>process.exit(1))"

# tini is PID 1 so orphaned tool processes (agents start them in sessions of their own) are reaped.
ENTRYPOINT ["/usr/bin/tini", "--", "/app/entrypoint.sh"]
