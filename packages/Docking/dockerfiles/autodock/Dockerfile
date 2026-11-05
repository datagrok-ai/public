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

# VEX: the runtime autodock env is Python 2.7 — its flask/werkzeug/jinja2/numpy have no
# py2.7 upgrades and are tracked in VEX. Its install-only pip/setuptools/wheel/certifi are
# removed, as are the py2.7 (and base) ensurepip bootstrap wheels; the floors of the
# build-only base (py3.14) env, where the leftover awscli deps live, are raised; its pip is
# removed too, as the vendored copies it ships (urllib3, msgpack, setuptools) are flagged.
RUN mamba remove -n autodock --force -y pip setuptools wheel certifi && \
    python -c "import flask, concurrent.futures; from MolKit import Read; from AutoDockTools.MoleculePreparation import AD4LigandPreparation, AD4ReceptorPreparation" && \
    rm -f /usr/local/miniforge/envs/autodock/lib/python2.7/ensurepip/_bundled/*.whl \
          /usr/local/miniforge/lib/python3.*/ensurepip/_bundled/*.whl && \
    /usr/local/miniforge/bin/pip install --no-cache-dir --upgrade 'urllib3>=2.8' 'msgpack>=1.2.3' && \
    mamba install -n base -c conda-forge -y 'py-rattler>=0.28.1' && \
    mamba remove -n base --force -y pip && \
    mamba clean -afy

COPY autodock.py /app/autodock.py

EXPOSE 8000

ENTRYPOINT ["python", "autodock.py", "--host", "0.0.0.0", "--port", "8000"]