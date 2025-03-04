FROM python:3.10
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
    curl \
    git \
    wget \
    unzip \
    bzip2 \
    xz-utils \
    libgomp1 \
    swig \
    cmake \
    build-essential \
    python3-dev \
    libeigen3-dev \
    zlib1g-dev \
    libffi-dev \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY . /app
RUN if [ "$(uname -m)" = "aarch64" ]; then \
        curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -o /tmp/miniforge.sh; \
    else \
        curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o /tmp/miniforge.sh; \
    fi && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh
USER root
ENV PATH="/app/venv/bin:$PATH"
RUN python3.10 -m venv /app/venv && \
    /app/venv/bin/pip install --upgrade pip && \
    /app/venv/bin/pip install aizynthfinder[all]
RUN mkdir ./aizynthcli_data
RUN /app/venv/bin/python -c "import sys; import aizynthfinder; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', 'aizynthcli_data']; main()"
RUN /app/venv/bin/pip install Flask
RUN /app/venv/bin/pip install flask_cors
EXPOSE 8000
ENTRYPOINT ["python", "app.py", "--host", "0.0.0.0", "--port", "8000"]