# Add ability to easily change python version
ARG PYTHON_VERSION=3.8

# Use our base image which is optimised and secured by our flow
FROM datagrok/python:${PYTHON_VERSION}

MAINTAINER Volodymyr Dyma <vdyma@datagrok.ai>

# Use pipefail to avoid errors after piping
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG HOME_DIR=/home/grok
# It is important to run docker container as user and not as root
RUN groupadd --gid 2001 grok ; \
    useradd --system --create-home --home ${HOME_DIR} --gid grok --uid 1001 grok

# Remove unzip after unpacking
ARG DEBIAN_FRONTEND=noninteractive
RUN savedAptMark="$(apt-mark showmanual)" ; \
    apt-get update ; \
    apt-get install -y --no-install-recommends \
      unzip \
      wget \
    ; \
    apt-mark auto '.*' > /dev/null ; \
    [ -z "$savedAptMark" ] || apt-mark manual $savedAptMark ; \
    wget https://mafft.cbrc.jp/alignment/software/mafft_7.520-1_amd64.deb -O mafft.deb; \
    apt install -y ./mafft.deb; \
    rm -rf mafft.deb; \
    wget https://github.com/Merck/PepSeA/archive/refs/heads/main.zip -O PepSeA.zip; \
    unzip -q PepSeA.zip -d /opt ; \
    rm -rf PepSeA.zip; \
    chown -R grok:grok /opt/PepSeA-main ; \
    apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false ; \
    apt-get clean ; \
    rm -rf /var/lib/apt/lists/* /var/cache/* /var/tmp/* /tmp/* /opt/PepSeA-main/alignment/mafft /opt/PepSeA-main/local.env; \
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

# It is important to run docker container as user and not as root
USER grok:grok


# Command source does not work for Docker, cause it will apply only to one layer
# The PATH works better for Docker
ARG VIRTUAL_ENV=${HOME_DIR}/venv
ENV PATH "${VIRTUAL_ENV}/bin:$PATH"
RUN mkdir -p ${VIRTUAL_ENV} ; \
    python -m venv ${VIRTUAL_ENV} ; \
    pip install --no-cache-dir --upgrade pip ; \
    pip install --no-cache-dir --upgrade wheel setuptools ; \
    pip install --timeout 3600 --no-cache-dir \
      --requirement /opt/PepSeA-main/requirements.txt ; \
    find ${VIRTUAL_ENV} -type d \( -path *tensorflow \) -prune -false -depth \( \
      \( -type d -a \( -name test -o -name tests -o -name idle_test \) \) \
          -o \( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
      \) -exec rm -rf '{}' + ; \
    rm -rf /var/cache/* /var/tmp/* /tmp/*

WORKDIR /opt/PepSeA-main

ENV ADD_HTTPS_MIDDLEWARE=""
ENV MAFFT_DIR="/usr/bin/"

# CMD can be changed easially, you need ENTRYPOINT in this case
ENTRYPOINT ["uvicorn", "--host", "0.0.0.0", "alignment.api:api"]

EXPOSE 8000
