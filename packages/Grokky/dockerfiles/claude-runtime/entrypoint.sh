#!/bin/sh
set -e

setpriv --reuid=1002 --regid=2002 --clear-groups env HOME=/home/broker node /app/dist/broker/broker.js &
BROKER_PID=$!

i=0
while [ "$i" -lt 50 ]; do
  if node -e "fetch('http://127.0.0.1:8377/health').then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))" 2>/dev/null; then
    break
  fi
  i=$((i + 1))
  sleep 0.2
done

exec setpriv --reuid=1001 --regid=2001 --clear-groups \
  env -i \
    PATH="$PATH" \
    HOME=/home/grok \
    CLAUDE_WORKSPACE="${CLAUDE_WORKSPACE:-/workspace}" \
    DATAGROK_API_URL="${DATAGROK_API_URL:-}" \
    NODE_ENV="${NODE_ENV:-}" \
    BROKER_PID="$BROKER_PID" \
    node /app/dist/server.js
