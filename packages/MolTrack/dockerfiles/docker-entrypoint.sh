#!/bin/bash
set -e

echo "Parsing DB environment variable..."
echo "$db"

if [ -n "$db" ]; then
  DB_USER=$(echo "$db" | jq -r '.credentials.parameters.login')
  DB_PASSWORD=$(echo "$db" | jq -r '.credentials.parameters.password')
  DB_NAME=$(echo "$db" | jq -r '.parameters.db')
  APP_PORT=28000

  DB_HOST_RAW=$(echo "$db" | jq -r '.parameters.server')
  echo $DB_HOST_RAW

  if [[ "$DB_HOST_RAW" == "localhost" ]]; then
    DB_HOST='host.docker.internal'
  else
    DB_HOST=$DB_HOST_RAW
  fi

  export DB_HOST DB_USER DB_PASSWORD DB_NAME APP_PORT

else
  echo "Environment variable 'db' not set."
  exit 1
fi

exec moltrack
