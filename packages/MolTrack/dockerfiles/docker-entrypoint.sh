#!/bin/bash
set -e

echo "Parsing DB environment variable..."

if [ -n "$db" ]; then
  export DB_HOST="host.docker.internal"
  export DB_USER=$(echo "$db" | jq -r '.credentials.parameters.login')
  export DB_PASSWORD=$(echo "$db" | jq -r '.credentials.parameters.password')
  export DB_NAME=$(echo "$db" | jq -r '.parameters.db')
  export APP_PORT=18000
else
  echo "Environment variable 'db' not set."
  exit 1
fi

exec moltrack