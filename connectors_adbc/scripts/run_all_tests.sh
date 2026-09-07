#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

echo "=== Stage 1: Unit Tests ==="
cargo test -- --test-threads=4

# Needs a local ClickHouse (docker container `clickhouse`, localhost:8123) and the
# adbc_clickhouse driver in lib/ or $ADBC_CLICKHOUSE_DRIVER.
echo "=== Stage 2: Integration Tests — ClickHouse ==="
cargo test --test clickhouse_integration -- --ignored --test-threads=1

echo "=== Stage 3: E2E Tests (in-process server) ==="
cargo test --test e2e_http -- --ignored --test-threads=1

echo "=== Stage 4: Docker Build ==="
docker build -t grok_connect_adbc:test .

# The server takes no credential config — the smoke tests send credentials in the request
# payload, reading them from DATABRICKS_PAT_PATH / DATABRICKS_JDBC_URL / SNOWFLAKE_CREDS_PATH
# on this host. An unset variable skips the corresponding test.
echo "=== Stage 5: E2E Smoke Tests (Docker container) ==="
docker network create grok_connect_adbc_test_net 2>/dev/null || true
docker run -d --name grok_connect_adbc_smoke \
  --network grok_connect_adbc_test_net \
  -p 1234:1234 \
  grok_connect_adbc:test
trap 'docker rm -f grok_connect_adbc_smoke 2>/dev/null; docker network rm grok_connect_adbc_test_net 2>/dev/null' EXIT
bash scripts/wait_for_service.sh http://localhost:1234/health 60
bash scripts/e2e_smoke.sh http://localhost:1234
docker rm -f grok_connect_adbc_smoke

echo "=== Stage 6: Platform Integration Test (optional) ==="
if [ "${SKIP_PLATFORM_TEST:-0}" = "1" ]; then
  echo "SKIPPED (SKIP_PLATFORM_TEST=1)"
else
  bash scripts/platform_test.sh || {
    echo "WARNING: Platform test failed — this is non-blocking for PoC"
    echo "Set SKIP_PLATFORM_TEST=1 to skip."
  }
fi

echo "=== ALL TESTS PASSED ==="
