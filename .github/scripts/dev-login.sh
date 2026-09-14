#!/usr/bin/env bash
# Exchanges a developer key for a session token and echoes it.
#
# Usage: dev-login.sh <apiUrl> <devKey>
#
# The key goes in the Authorization header. Servers before 1.28 only knew the key-in-URL
# form and answer 401 to the key-less route (it is not on their anonymous allow-list), so
# the old form is tried next — a workflow run can meet an older datagrok image than the
# checkout it came from.
set -euo pipefail

apiUrl="${1%/}"
key="$2"

token="$(curl -s -X POST "${apiUrl}/users/login/dev" -H "Authorization: Dev ${key}" | jq -r '.token // empty')"
if [ -z "$token" ]; then
  token="$(curl -s -X POST "${apiUrl}/users/login/dev/${key}" | jq -r '.token // empty')"
fi

echo "$token"
