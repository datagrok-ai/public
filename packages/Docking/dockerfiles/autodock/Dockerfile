FROM cgr.dev/chainguard/wolfi-base AS build
RUN apk add --no-cache build-base opencl-dev curl unzip
WORKDIR /build

# AutoDockTools 1.5.7 ported to Python 3 (MolKit, PyBabel, mglutil): prepare_ligand4 / prepare_receptor4
ARG ADT_PY3_SHA=a62a5d98116f1590183b58a9ad732b997cf2579c
RUN curl -fsSL https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3/archive/${ADT_PY3_SHA}.tar.gz | tar -xz \
    && mkdir -p /out/adt \
    && for p in AutoDockTools MolKit PyBabel mglutil; do cp -r AutoDockTools_py3-${ADT_PY3_SHA}/$p /out/adt/; done

ARG CLINFO_TAG=3.0.25.02.14
RUN curl -fsSL https://github.com/Oblomov/clinfo/archive/refs/tags/${CLINFO_TAG}.tar.gz | tar -xz \
    && make -C clinfo-${CLINFO_TAG} \
    && install -D -m 755 clinfo-${CLINFO_TAG}/clinfo /out/root/usr/local/bin/clinfo

RUN curl -fsSL https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-x86_64Linux2.tar -o autodocksuite.tar \
    && echo "d9f8bada0b0618da1b2c9e58aa65188bbde3ef5c36d278b7b8657169502cbb03  autodocksuite.tar" | sha256sum -c - \
    && tar -xf autodocksuite.tar x86_64Linux2/autogrid4 \
    && install -D -m 755 x86_64Linux2/autogrid4 /out/root/usr/local/x86_64Linux2/autogrid4 \
    && curl -fsSL https://github.com/ccsb-scripps/AutoDock-GPU/releases/download/v1.5.3/adgpu-v1.5.3_linux_ocl_128wi -o autodock-gpu \
    && echo "05a80fd91016d8c4e7b9e9072bc3a2071a83e16427066329f4e88a88a2fee25d  autodock-gpu" | sha256sum -c - \
    && install -D -m 755 autodock-gpu /out/root/opt/autodock-gpu

# Precomputed grids of the demo targets
RUN s3=https://datagrok-data.s3.us-east-2.amazonaws.com && mkdir -p /out/root/app \
    && for key in $(curl -fsSL "$s3/?list-type=2&prefix=autodock/" | grep -o '<Key>[^<]*\.zip</Key>' | sed -e 's/<Key>//' -e 's/<\/Key>//'); do \
         curl -fsSL "$s3/$key" -o grids.zip && unzip -oq grids.zip -x '__MACOSX/*' '*.DS_Store' -d /out/root/app && rm grids.zip; \
       done \
    && ls /out/root/app

FROM cgr.dev/chainguard/wolfi-base
RUN apk upgrade --no-cache \
    && apk add --no-cache python-3.13 py3.13-flask py3.13-numpy-2.4 libstdc++ libgcc libgomp opencl \
    && rm -f /usr/bin/clinfo

RUN mkdir -p /etc/OpenCL/vendors && \
    echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

COPY --from=build /out/adt/ /usr/lib/python3.13/site-packages/
COPY --from=build /out/root/ /
RUN python3 -m compileall -q /usr/lib/python3.13/site-packages/AutoDockTools /usr/lib/python3.13/site-packages/MolKit \
      /usr/lib/python3.13/site-packages/PyBabel /usr/lib/python3.13/site-packages/mglutil > /dev/null 2>&1; \
    for b in /opt/autodock-gpu /usr/local/x86_64Linux2/autogrid4 /usr/local/bin/clinfo; do \
      if /lib64/ld-linux-x86-64.so.2 --list $b | grep 'not found'; then exit 1; fi; \
    done

WORKDIR /app
COPY autodock.py /app/autodock.py

RUN --mount=type=bind,source=selftest,target=/selftest python3 /selftest/selftest.py

EXPOSE 8000

ENTRYPOINT ["python3", "autodock.py", "--host", "0.0.0.0", "--port", "8000"]
