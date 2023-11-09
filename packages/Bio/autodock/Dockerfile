FROM python:2.7

ENV PYTHONUNBUFFERED 1

WORKDIR /app

RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y ocl-icd-libopencl1

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh \
    && /bin/bash /tmp/miniforge.sh -b -p /usr/local/miniforge \
    && rm /tmp/miniforge.sh

ENV PATH="/usr/local/miniforge/bin:$PATH"

RUN conda install -c conda-forge mamba

RUN mamba create -n autodock python=2.7
RUN echo "source activate autodock" > ~/.bashrc
ENV PATH /usr/local/miniforge/envs/autodock/bin:$PATH
RUN mamba install -n autodock -c insilichem autodocktools-prepare

RUN wget https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-x86_64Linux2.tar -O /tmp/autodocksuite.tar \
    && tar -xvf /tmp/autodocksuite.tar -C /usr/local/

RUN wget https://github.com/ccsb-scripps/AutoDock-GPU/releases/download/v1.5.3/adgpu-v1.5.3_linux_ocl_128wi -O /opt/autodock-gpu
RUN chmod +x /opt/autodock-gpu

RUN pip install Flask

COPY autodock.py /app/autodock.py

EXPOSE 8000

ENTRYPOINT ["python", "autodock.py", "--host", "0.0.0.0", "--port", "8000"]