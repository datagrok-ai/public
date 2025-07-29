FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    DB_USER=postgres \
    DB_PASSWORD=postgres \
    DB_HOST=172.17.0.1 \
    DB_PORT=5432 \
    DB_NAME=moltrack \
    DB_SCHEMA=moltrack

RUN pip install --no-cache-dir dg-mol-track==0.2.0

EXPOSE 8000

ENTRYPOINT ["moltrack"]