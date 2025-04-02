FROM python:3.10

ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies
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
    libgomp1 \
    libeigen3-dev \
    zlib1g-dev \
    libffi-dev \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up application
WORKDIR /app
ENV PATH="/app/venv/bin:$PATH"

# Create virtual environment and install packages
RUN python3.10 -m venv /app/venv && \
    /app/venv/bin/pip install --upgrade pip && \
    /app/venv/bin/pip install aizynthfinder[all] Flask flask_cors gunicorn

# Download public data
RUN mkdir ./aizynthcli_data && \
    /app/venv/bin/python -c "import sys; import aizynthfinder; from aizynthfinder.tools.download_public_data import main; sys.argv = ['download_public_data', 'aizynthcli_data']; main()"

# Copy application
COPY app.py /app/

EXPOSE 8000

# Run with Gunicorn (--preload is critical for multiprocessing.Manager)
CMD ["gunicorn", "--preload", "-w", "4", "-b", "0.0.0.0:8000", "--timeout", "300", "app:app"]