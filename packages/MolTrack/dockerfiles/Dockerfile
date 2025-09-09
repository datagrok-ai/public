FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y jq && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir dg-mol-track==0.6.11

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh

RUN sed -i 's/\r$//' /usr/local/bin/docker-entrypoint.sh

RUN chmod +x /usr/local/bin/docker-entrypoint.sh

EXPOSE 18000

ENTRYPOINT ["docker-entrypoint.sh"]