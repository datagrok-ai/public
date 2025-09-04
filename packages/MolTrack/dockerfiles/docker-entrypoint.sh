#!/bin/bash
set -e

echo "Parsing DB environment variable..."

echo $db

if [ -n "$db" ]; then
  # export DB_HOST=$(echo "$db" | jq -r '.parameters.server')
  export DB_HOST="10.152.183.177"
  export DB_USER=$(echo "$db" | jq -r '.credentials.parameters.login')
  export DB_PASSWORD=$(echo "$db" | jq -r '.credentials.parameters.password')
  export DB_NAME=$(echo "$db" | jq -r '.parameters.db')
  export APP_PORT=18000
else
  echo "Environment variable 'db' not set."
  exit 1
fi

exec moltrack