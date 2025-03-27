FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y \
    apt-utils \
    unzip \
    tar \
    curl \
    xz-utils \
    ocl-icd-libopencl1 \
    opencl-headers \
    clinfo \
    libgomp1 \
    wget \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility

RUN apt-get update && apt-get install -y clinfo

WORKDIR /app

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

RUN curl -LJO https://raw.githubusercontent.com/MolecularFlipbook/FlipbookApp/master/mfb/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py \
    && mv prepare_dpf4.py /usr/local/bin/ \
    && chmod +x /usr/local/bin/prepare_dpf4.py

RUN pip install Flask

RUN pip install futures

RUN mamba install -n base -c conda-forge awscli -y && \
    mamba run -n base aws s3 cp s3://datagrok-data/autodock /app --no-sign-request --recursive --region us-east-2 && \
    find /app -name '*.zip' -exec unzip -o {} -d /app \; && \
    rm -rf /app/*.zip && \
    mamba remove -n base awscli -y && \
    mamba clean -afy

COPY autodock.py /app/autodock.py

EXPOSE 8000

ENTRYPOINT ["python", "autodock.py", "--host", "0.0.0.0", "--port", "8000"]