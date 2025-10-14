FROM python:3.11-slim

SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y --no-install-recommends \
    curl bzip2 ca-certificates gcc git \
    && apt clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

RUN python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('all-MiniLM-L6-v2')"

COPY . /app

EXPOSE 8000

ENTRYPOINT ["python", "-m", "app", "--host", "0.0.0.0", "--port", "8000"]