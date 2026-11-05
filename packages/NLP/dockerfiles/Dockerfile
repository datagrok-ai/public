# Wolfi (glibc, near-zero CVEs): torch ships manylinux wheels only, so no musl base.
FROM cgr.dev/chainguard/wolfi-base AS build

RUN apk add --no-cache python-3.11 py3.11-pip

COPY requirements.txt /tmp/requirements.txt
RUN python3.11 -m venv --without-pip /opt/venv && \
    python3.11 -m pip install --no-cache-dir --python /opt/venv/bin/python -r /tmp/requirements.txt

ENV PATH=/opt/venv/bin:$PATH
RUN python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('all-MiniLM-L6-v2')"

FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache && apk add --no-cache python-3.11 libstdc++ libgomp && \
    rm -rf /usr/lib/python3.11/ensurepip

ENV PATH=/opt/venv/bin:$PATH
WORKDIR /app

COPY --from=build /opt/venv /opt/venv
COPY --from=build /root/.cache/huggingface /root/.cache/huggingface
COPY . /app

RUN HF_HUB_OFFLINE=1 python -c "import app; from sentence_transformers import SentenceTransformer; \
assert len(SentenceTransformer('all-MiniLM-L6-v2').encode('hello')) == 384" && rm -f /app/app.log

EXPOSE 8000

ENTRYPOINT ["python", "-m", "app", "--host", "0.0.0.0", "--port", "8000"]
