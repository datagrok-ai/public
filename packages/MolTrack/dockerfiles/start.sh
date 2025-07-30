#!/bin/bash
set -e

DATA_DIR=/var/lib/postgresql/data
POSTGRES_BIN=/usr/lib/postgresql/15/bin

if [ ! -s "$DATA_DIR/PG_VERSION" ]; then
  echo "Initializing DB cluster..."
  gosu postgres $POSTGRES_BIN/initdb -D "$DATA_DIR"
fi

echo "Starting Postgres..."
gosu postgres $POSTGRES_BIN/postgres -D "$DATA_DIR" &

echo "Waiting for Postgres..."
until pg_isready -h localhost -p 5432; do
  sleep 1
done

echo "Postgres ready."

echo "Running 01_db.sql on 'postgres' database..."
gosu postgres psql -v ON_ERROR_STOP=1 --username postgres --dbname postgres -f /docker-entrypoint-initdb.d/01_db.sql

echo "Running other SQL scripts on 'moltrack' database..."
for sqlfile in /docker-entrypoint-initdb.d/*.sql; do
  case "$sqlfile" in
    */01_db.sql) continue ;;
    *)
      echo "Running $sqlfile"
      gosu postgres psql -v ON_ERROR_STOP=1 --username postgres --dbname moltrack -f "$sqlfile"
      ;;
  esac
done

touch "$DATA_DIR/.db_init_done"

echo "Starting moltrack..."
/opt/venv/bin/moltrack