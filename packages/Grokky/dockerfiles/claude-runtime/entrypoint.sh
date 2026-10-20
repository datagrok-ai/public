#!/bin/sh
set -e

# Started as anyone but root, the server could not drop to uid 1001 itself and would stay dumpable.
[ "$(id -u)" = 0 ] || { echo "Start the Grokky runtime as root: it drops to uid 1001 itself." >&2; exit 1; }

# Only this allowlist of the container's environment reaches the server, so no foreign or stale credential does.
exec env -i \
    PATH="$PATH" \
    HOME=/home/grok \
    CLAUDE_WORKSPACE="${CLAUDE_WORKSPACE:-/workspace}" \
    DATAGROK_API_URL="${DATAGROK_API_URL:-}" \
    provider="${provider:-}" \
    apiKey="${apiKey:-}" \
    baseUrl="${baseUrl:-}" \
    region="${region:-}" \
    apiFormat="${apiFormat:-}" \
    model="${model:-}" \
    fastModel="${fastModel:-}" \
    TEST_UPSTREAM="${TEST_UPSTREAM:-}" \
    HERMES_MAX_PROCESSES="${HERMES_MAX_PROCESSES:-}" \
    HERMES_IDLE_TTL_MS="${HERMES_IDLE_TTL_MS:-}" \
    HERMES_WARM_SPARES="${HERMES_WARM_SPARES:-}" \
    HERMES_READY_TIMEOUT_MS="${HERMES_READY_TIMEOUT_MS:-}" \
    HERMES_MAX_SESSIONS="${HERMES_MAX_SESSIONS:-}" \
    HERMES_COMPRESSION_THRESHOLD_TOKENS="${HERMES_COMPRESSION_THRESHOLD_TOKENS:-}" \
    HERMES_LOG_STDERR="${HERMES_LOG_STDERR:-}" \
    node /app/dist/server.js
