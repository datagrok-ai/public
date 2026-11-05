ARG PYTHON_VERSION=3.11
ARG HOME_DIR=/home/grok
ARG VIRTUAL_ENV=${HOME_DIR}/venv
ARG MAFFT_VERSION=7.525

# Wolfi (glibc, near-zero CVEs): mafft is built from source and the venv is built here,
# then only they, PepSeA and the interpreter reach the final stage
FROM cgr.dev/chainguard/wolfi-base AS builder

ARG PYTHON_VERSION
ARG VIRTUAL_ENV
ARG MAFFT_VERSION

RUN apk add --no-cache bash build-base curl unzip python-${PYTHON_VERSION} py${PYTHON_VERSION}-pip

SHELL ["/bin/bash", "-e", "-o", "pipefail", "-c"]

RUN curl -fsSL https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz -o /tmp/mafft.tgz ; \
    echo "edb34ae9b26d6b55328c18fa060ed741bca8cd599c2f4f8fad0e0871c8082265  /tmp/mafft.tgz" | sha256sum -c - ; \
    tar -xzf /tmp/mafft.tgz -C /tmp ; \
    # mafft's K&R-style prototypes do not compile as C23, gcc's default since 15 \
    make -C /tmp/mafft-${MAFFT_VERSION}-without-extensions/core PREFIX=/usr/lib/mafft CFLAGS='-O3 -std=gnu17 -fpermissive' install ; \
    rm -rf /usr/lib/mafft/share /usr/lib/mafft/bin/*.rb /tmp/*

# PepSeA's requirements.txt pins 2021 versions that do not install on python 3.11; only what alignment/ imports is installed.
# pip stays in the builder: the venv is created without it
RUN python${PYTHON_VERSION} -m venv --without-pip ${VIRTUAL_ENV} ; \
    python${PYTHON_VERSION} -m pip --python ${VIRTUAL_ENV}/bin/python install --timeout 3600 --no-cache-dir \
      'fastapi>=0.133' 'starlette>=1.3.1' 'pydantic>=2.7,<3' \
      'uvicorn[standard]>=0.34' 'h11>=0.16' 'ujson>=5.12.1' 'pymsa==0.8.1' ; \
    find ${VIRTUAL_ENV} -depth \( \
      \( -type d -a \( -name test -o -name tests -o -name idle_test \) \) \
          -o \( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
      \) -exec rm -rf '{}' +

RUN curl -fsSL https://github.com/Merck/PepSeA/archive/refs/heads/main.zip -o /tmp/PepSeA.zip ; \
    unzip -q /tmp/PepSeA.zip -d /opt ; \
    rm -rf /tmp/PepSeA.zip /opt/PepSeA-main/alignment/mafft /opt/PepSeA-main/local.env ; \
    sed -i "s/files = glob.glob(f'\*{timestamp}\*')/files = glob.glob(f'\*{timestamp}\*.txt')/" /opt/PepSeA-main/alignment/AlignSubPeptides.py; \
    echo $'\n\
@api.get("/distout", response_class=PlainTextResponse)\n\
def distout():\n\
    """Returns the output of --distout parameter as a plain text."""\n\
    \n\
    distout_output = None\n\
\n\
    for file in os.listdir("./"):\n\
        if not file.endswith(".hat2"):\n\
            continue\n\
\n\
        path = os.path.join("./", file)\n\
        with open(path, "r", encoding="utf-8") as f:\n\
            distout_output = "".join(f.readlines())\n\
        os.remove(path)\n\
        break\n\
\n\
    if distout_output is None:\n\
        raise FileNotFoundError("distance matrix is not found")\n\
\n\
    return distout_output\n\
\n\
from fastapi import Request, HTTPException\n\
import sys\n\
import ujson as json\n\
from io import StringIO\n\
from subprocess import CalledProcessError\n\
\n\
@api.middleware("http")\n\
async def error_to_json_middleware(request: Request, call_next):\n\
    sys.stderr = temp_stderr = StringIO()\n\
    response = None\n\
    try:\n\
        response = await call_next(request)\n\
    finally:\n\
        sys.stderr = sys.__stderr__\n\
        stderr_msg = temp_stderr.getvalue()\n\
        if len(stderr_msg) > 0 and response and response.headers["content-type"] == "application/json":\n\
            response_body = [chunk async for chunk in response.body_iterator]\n\
            response_str = "".join(((s.decode() if isinstance(s, bytes) else str(s)) for s in response_body))\n\
            response_obj = json.loads(response_str)\n\
            response_obj.update({"pepsea-stderr": stderr_msg})\n\
            headers = {}\n\
            for header_name in response.headers:\n\
                if header_name != "content-length" and header_name != "content-type":\n\
                    headers[header_name] = response.headers[header_name]\n\
            json_response = JSONResponse(status_code=response.status_code, headers=headers, content=response_obj)\n\
            return json_response\n\
    return response\n\
\n\
\n\
def pepsea_json_response(err_msg: str):\n\
    return JSONResponse(status_code=500, content={\n\
        "pepsea-error": err_msg,\n\
    })\n\
\n\
\n\
@api.exception_handler(CalledProcessError)\n\
async def called_process_error_handler(request: Request, ex: CalledProcessError):\n\
    s = ex.stderr\n\
    err_msg = s.decode() if isinstance(s, bytes) else str(ex.stderr)\n\
    sys.__stderr__.write("called_process_error_handler() {0}: {1}\\n".format(type(ex).__name__, err_msg))\n\
    return pepsea_json_response(err_msg)\n\
\n\
\n\
@api.exception_handler(Exception)\n\
async def base_exception_handler(request: Request, ex: Exception):\n\
    err_msg = str(ex)\n\
    sys.__stderr__.write("base_exception_handler() {0}: {1}\\n".format(type(ex).__name__, err_msg))\n\
    return pepsea_json_response(err_msg)\n\
\n' >> /opt/PepSeA-main/alignment/api.py;

FROM cgr.dev/chainguard/wolfi-base

LABEL maintainer="Volodymyr Dyma <vdyma@datagrok.ai>"

ARG PYTHON_VERSION
ARG HOME_DIR
ARG VIRTUAL_ENV

# bash: the mafft driver scripts are bash
RUN set -e ; \
    apk upgrade --no-cache ; \
    apk add --no-cache bash python-${PYTHON_VERSION} ; \
    rm -rf /usr/lib/python3*/ensurepip ; \
    addgroup -g 2001 grok ; \
    adduser -D -h ${HOME_DIR} -G grok -u 1001 grok

COPY --from=builder /usr/lib/mafft /usr/lib/mafft
COPY --from=builder --chown=grok:grok ${VIRTUAL_ENV} ${VIRTUAL_ENV}
COPY --from=builder --chown=grok:grok /opt/PepSeA-main /opt/PepSeA-main

# It is important to run docker container as user and not as root
USER grok:grok

ENV PATH="${VIRTUAL_ENV}/bin:$PATH"

WORKDIR /opt/PepSeA-main

ENV ADD_HTTPS_MIDDLEWARE=""
ENV MAFFT_DIR="/usr/lib/mafft/bin/"

RUN "${MAFFT_DIR}ginsi" --version && python -c "import alignment.api, uvloop, httptools, ujson, pymsa" && \
    printf '>a\nMKVLA\n>b\nMKVA\n' > /tmp/t.fa && "${MAFFT_DIR}mafft" --quiet /tmp/t.fa && rm -f /tmp/t.fa

# CMD can be changed easially, you need ENTRYPOINT in this case
ENTRYPOINT ["uvicorn", "--host", "0.0.0.0", "alignment.api:api"]

EXPOSE 8000
