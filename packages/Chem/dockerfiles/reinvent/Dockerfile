FROM python:3.10

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
    curl \
    wget \
    git \
    unzip \
    tar \
    bzip2 \
    xz-utils \
    libgomp1 \
    swig \
    cmake \
    build-essential \
    python3-dev \
    libeigen3-dev \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN git clone https://github.com/MolecularAI/REINVENT4.git && \
    git clone https://github.com/MolecularAI/DockStream.git

RUN find /app/REINVENT4 -type f -exec sed -i 's/model_dump()/model_config/g' {} +

COPY . /app

RUN python3.10 -m venv /app/venv && \
    /app/venv/bin/pip install --upgrade pip && \
    /app/venv/bin/pip install -r requirements-macOS.lock && \
    cd /app/REINVENT4 && \
    /app/venv/bin/pip install --no-deps .

RUN if [ "$(uname -m)" = "aarch64" ]; then \
        curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -o /tmp/miniforge.sh; \
    else \
        curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o /tmp/miniforge.sh; \
    fi && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

USER root

ENV PATH="/opt/conda/bin:$PATH"

RUN conda update -n base -c defaults conda -y

WORKDIR /app/DockStream

RUN sed -i 's/python = 3.7/python = 3.10/' environment.yml

RUN conda env create -f environment.yml && \
    conda clean --all --yes

RUN chmod -R 777 /opt/conda/envs/DockStream

WORKDIR /app

RUN mkdir -p /app/autodock_vina && \
    curl -L https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz -o /app/autodock_vina.tar.gz && \
    tar -xzvf /app/autodock_vina.tar.gz -C /app/autodock_vina && \
    rm /app/autodock_vina.tar.gz

RUN /app/venv/bin/pip install Flask

RUN /app/venv/bin/pip install descriptastorus

ENV PATH="/app/venv/bin:$PATH"

EXPOSE 6666

ENTRYPOINT ["python", "reinvent.py", "--host", "0.0.0.0", "--port", "6666"]