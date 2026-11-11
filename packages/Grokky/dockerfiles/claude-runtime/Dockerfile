ARG HERMES_COMMIT=1f6c02c07c697026888f1c66e634534fa2aae358
ARG NODE_IMAGE=node:22.23.4-trixie-slim@sha256:8eb619e402eaede1748a3c8173bffab391592812e85136452f781741c7242ddf

# Hermes' own package manager builds the venv, as Hermes' Dockerfile does: CPython and uv come from
# pm/lock.json (sha256-verified) into HERMES_RUNTIME_DIR, which the venv links to, so the path is kept.
# uv is needed for the build only: it first lifts the venv's libraries Hermes pins below their CVE fixes.
# Of the unstripped CPython the lock pins, only libpython is stripped: GNU strip
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
      --out /opt/hermes/.venv --sealed --extra anthropic --extra mcp --extra bedrock && \
    echo /opt/hermes > "$(/opt/hermes/.venv/bin/python -c 'import site; print(site.getsitepackages()[0])')/hermes_project.pth" && \
    uv="$(find "$HERMES_RUNTIME_DIR" -type f -name uv | head -n 1)" && [ -x "$uv" ] && \
    "$uv" pip install --python /opt/hermes/.venv/bin/python --upgrade \
      "pillow>=13.0.0" "pyjwt>=2.14.0" "certifi>=2026.10.18" "urllib3>=2.8.0" "mcp>=2.2.0" \
      $("$uv" pip list --python /opt/hermes/.venv/bin/python --format freeze | grep -io '^opentelemetry-[a-z-]*') && \
    rm -rf "$HERMES_RUNTIME_DIR"/uv-* \
      "$HERMES_RUNTIME_DIR"/python-*/lib/python3*/ensurepip \
      "$HERMES_RUNTIME_DIR"/python-*/lib/python3*/site-packages/pip \
      "$HERMES_RUNTIME_DIR"/python-*/lib/python3*/site-packages/pip-* && \
    strip --strip-debug "$HERMES_RUNTIME_DIR"/python-*/lib/libpython3.*.so.1.0 && \
    git rev-parse HEAD > COMMIT && \
    rm -rf .git tests website apps web ui-tui optional-skills && \
    /opt/hermes/.venv/bin/python -c "import run_agent, anthropic, mcp, jsonschema, ruamel.yaml, PIL.Image, jwt"

FROM ${NODE_IMAGE} AS app
WORKDIR /app
COPY package.json package-lock.json* tsconfig.json ./
COPY src ./src
RUN npm ci && npx tsc && rm -rf src tsconfig.json && npm prune --omit=dev

# Wolfi (glibc) carries no unfixed distro CVEs. The userland is what the agent's terminal and `!cmd` use;
# node is the official build from the Node image. npm and virtualenv are not shipped: nothing here runs them.
FROM cgr.dev/chainguard/wolfi-base
ARG HERMES_COMMIT
LABEL ai.datagrok.hermes.commit=$HERMES_COMMIT

RUN apk upgrade --no-cache && \
    apk add --no-cache bash coreutils findutils grep sed diffutils gzip git ca-certificates-bundle \
      python-3.13 py3.13-pip ripgrep tini libstdc++ libgcc && \
    test -x /usr/bin/tini
COPY --from=app /usr/local/bin/node /usr/local/bin/node

WORKDIR /app

ARG HOME_DIR=/home/grok
RUN addgroup -g 2001 grok && \
    adduser -D -h ${HOME_DIR} -G grok -u 1001 -s /bin/bash grok && \
    git clone --depth 1 --filter=blob:none --sparse https://github.com/datagrok-ai/public.git /workspace && \
    git -C /workspace sparse-checkout set --no-cone '/*' '!/docusaurus/' \
      '!*.gif' '!*.png' '!*.PNG' '!*.jpg' '!*.jpeg' '!*.bmp' '!*.mp4' '!*.jar' && \
    git clone --depth 1 --filter=blob:none --sparse https://github.com/datagrok-ai/help.git /workspace/help && \
    git -C /workspace/help sparse-checkout set --no-cone '/*' '!/uploads/' \
      '!*.gif' '!*.png' '!*.PNG' '!*.jpg' '!*.jpeg' '!*.bmp' '!*.mp4' && \
    mkdir -p /users /var/lib/grokky && \
    chown -R grok:grok ${HOME_DIR} /workspace /users /var/lib/grokky

COPY --from=app /app /app
COPY --from=hermes /opt/hermes /opt/hermes

COPY plugin /app/plugin
COPY agent /app/agent
COPY entrypoint.sh /app/entrypoint.sh
RUN sed -i 's/\r$//' /app/entrypoint.sh && chmod 0755 /app/entrypoint.sh && \
    node --version && python3 --version && pip3 --version && rg --version && git --version && \
    PYTHONPATH=/app/agent /opt/hermes/.venv/bin/python -c \
      "import run_agent, anthropic, mcp, jsonschema, ruamel.yaml, PIL.Image, jwt, grokky_agent.session"

ENV CLAUDE_WORKSPACE=/workspace

EXPOSE 5355

HEALTHCHECK --interval=5s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:5355/health').then(r=>r.ok?process.exit(0):process.exit(1)).catch(()=>process.exit(1))"

# tini is PID 1 so orphaned tool processes (agents start them in sessions of their own) are reaped.
ENTRYPOINT ["/usr/bin/tini", "--", "/app/entrypoint.sh"]
