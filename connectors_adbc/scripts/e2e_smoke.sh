#!/bin/bash
# E2E smoke tests using curl and jq against a running grok_connect_adbc instance.
# Usage: e2e_smoke.sh <base_url>
# Requires: curl, jq

set -euo pipefail

BASE_URL="${1:?Usage: e2e_smoke.sh <base_url>}"
FAILED=0

pass() { echo "  PASS: $1"; }
fail() { echo "  FAIL: $1"; FAILED=$((FAILED + 1)); }

echo "=== E2E Smoke Tests against $BASE_URL ==="

# Test 1: GET /health
echo "Test 1: GET /health"
RESP=$(curl -sf "$BASE_URL/health" 2>&1) && [ "$RESP" = "OK" ] \
  && pass "/health returns OK" \
  || fail "/health: expected OK, got: $RESP"

# Test 2: GET /info
echo "Test 2: GET /info"
NAME=$(curl -sf "$BASE_URL/info" | jq -r '.name' 2>&1)
[ "$NAME" = "grok_connect_adbc" ] \
  && pass "/info name is grok_connect_adbc" \
  || fail "/info: expected name=grok_connect_adbc, got: $NAME"

# Test 3: GET /conn
echo "Test 3: GET /conn"
CONN_RESP=$(curl -sf "$BASE_URL/conn")
echo "$CONN_RESP" | jq -e '.[] | select(.type == "Databricks")' > /dev/null 2>&1 \
  && pass "/conn lists Databricks" \
  || fail "/conn: Databricks not found"
echo "$CONN_RESP" | jq -e '.[] | select(.type == "Snowflake")' > /dev/null 2>&1 \
  && pass "/conn lists Snowflake" \
  || fail "/conn: Snowflake not found"
# Check canBrowseSchema=false for all
BROWSABLE=$(echo "$CONN_RESP" | jq '[.[] | .canBrowseSchema] | any' 2>&1)
[ "$BROWSABLE" = "false" ] \
  && pass "/conn all canBrowseSchema=false" \
  || fail "/conn: expected all canBrowseSchema=false"

# Test 4: POST /test (Databricks)
echo "Test 4: POST /test (Databricks)"
# DATABRICKS_PAT_PATH / DATABRICKS_JDBC_URL are required for this test; skipped if unset.
if [ -n "${DATABRICKS_PAT_PATH:-}" ] && [ -f "$DATABRICKS_PAT_PATH" ] && [ -n "${DATABRICKS_JDBC_URL:-}" ]; then
  PAT=$(cat "$DATABRICKS_PAT_PATH" | tr -d '[:space:]')
  JDBC_URL="$DATABRICKS_JDBC_URL"

  # Extract host and httpPath from JDBC URL
  HOST=$(echo "$JDBC_URL" | sed -n 's|.*://\([^:/]*\).*|\1|p')
  HTTP_PATH=$(echo "$JDBC_URL" | sed -n 's|.*httpPath=\([^;]*\).*|\1|p')

  DATABRICKS_CONN=$(cat <<EOJSON
{
  "dataSource": "Databricks",
  "parameters": {
    "server": "$HOST",
    "port": 443,
    "httpPath": "$HTTP_PATH",
    "db": "default"
  },
  "credentials": {
    "parameters": {
      "password": "$PAT"
    }
  }
}
EOJSON
)
  TEST_RESP=$(curl -sf -X POST -H "Content-Type: application/json" -d "$DATABRICKS_CONN" "$BASE_URL/test" 2>&1)
  echo "$TEST_RESP" | grep -q "Connection available" \
    && pass "Databricks connection test" \
    || fail "Databricks connection test: $TEST_RESP"
else
  echo "  SKIP: set DATABRICKS_PAT_PATH and DATABRICKS_JDBC_URL to run this test"
fi

# Test 5: POST /test (Snowflake)
echo "Test 5: POST /test (Snowflake)"
SF_CREDS_PATH="${SNOWFLAKE_CREDS_PATH:-}"
if [ -n "$SF_CREDS_PATH" ] && [ -f "$SF_CREDS_PATH" ]; then
  SF_SERVER=$(jq -r '.server // .Server' "$SF_CREDS_PATH")
  SF_DB=$(jq -r '.db // .Db // ""' "$SF_CREDS_PATH")
  SF_WH=$(jq -r '.warehouse // .Warehouse // ""' "$SF_CREDS_PATH")
  SF_LOGIN=$(jq -r '.login // .Login' "$SF_CREDS_PATH")
  SF_PASS=$(jq -r '.password // .Password' "$SF_CREDS_PATH")

  SNOWFLAKE_CONN=$(cat <<EOJSON
{
  "dataSource": "Snowflake",
  "parameters": {
    "server": "$SF_SERVER",
    "db": "$SF_DB",
    "warehouse": "$SF_WH"
  },
  "credentials": {
    "parameters": {
      "login": "$SF_LOGIN",
      "password": "$SF_PASS"
    }
  }
}
EOJSON
)
  TEST_RESP=$(curl -sf -X POST -H "Content-Type: application/json" -d "$SNOWFLAKE_CONN" "$BASE_URL/test" 2>&1)
  echo "$TEST_RESP" | grep -q "Connection available" \
    && pass "Snowflake connection test" \
    || fail "Snowflake connection test: $TEST_RESP"
else
  echo "  SKIP: set SNOWFLAKE_CREDS_PATH to run this test"
fi

echo ""
if [ "$FAILED" -gt 0 ]; then
  echo "=== $FAILED TEST(S) FAILED ==="
  exit 1
fi
echo "=== ALL SMOKE TESTS PASSED ==="
