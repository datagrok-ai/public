#!/bin/bash
# Poll a health endpoint until it returns 200 OK, with timeout.
# Usage: wait_for_service.sh <url> [timeout_seconds]

set -euo pipefail

URL="${1:?Usage: wait_for_service.sh <url> [timeout_seconds]}"
TIMEOUT="${2:-60}"

echo "Waiting for $URL (timeout: ${TIMEOUT}s)..."

ELAPSED=0
INTERVAL=2

while [ "$ELAPSED" -lt "$TIMEOUT" ]; do
  if curl -sf "$URL" > /dev/null 2>&1; then
    echo "Service ready at $URL (after ${ELAPSED}s)"
    exit 0
  fi
  sleep "$INTERVAL"
  ELAPSED=$((ELAPSED + INTERVAL))
done

echo "ERROR: Service at $URL not ready after ${TIMEOUT}s"
exit 1
