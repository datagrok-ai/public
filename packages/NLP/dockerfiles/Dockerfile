FROM python:3.11-slim

SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    rm -f /usr/local/lib/python3.11/ensurepip/_bundled/*.whl

RUN python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('all-MiniLM-L6-v2')"

RUN dpkg --purge --force-remove-essential --force-depends \
      perl perl-base perl-modules-5.40 libperl5.40

COPY . /app

EXPOSE 8000

ENTRYPOINT ["python", "-m", "app", "--host", "0.0.0.0", "--port", "8000"]