FROM cgr.dev/chainguard/wolfi-base AS pip_builder

ARG VIRTUAL_ENV=/app/venv

RUN apk add --no-cache python-3.12 py3.12-pip libstdc++

COPY requirements.in ${VIRTUAL_ENV}/requirements.in

RUN python3.12 -m venv --without-pip ${VIRTUAL_ENV} && \
    python3.12 -m pip --python ${VIRTUAL_ENV}/bin/python install --timeout 3600 --no-cache-dir \
      -r ${VIRTUAL_ENV}/requirements.in && \
    find ${VIRTUAL_ENV} -depth -type d \( -name __pycache__ -o -name tests \) -path '*/site-packages/*' -exec rm -rf {} +


# Wolfi (glibc) with only the interpreter: the venv is copied in without pip/setuptools/wheel.
FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache && \
    apk add --no-cache python-3.12 libstdc++ && \
    rm -rf /usr/lib/python3.12/ensurepip

ARG CHEM_USER=chem
RUN addgroup -S $CHEM_USER && adduser -S -G $CHEM_USER $CHEM_USER
USER $CHEM_USER

WORKDIR /app

ARG VIRTUAL_ENV=/app/venv

COPY --from=pip_builder --chown=$CHEM_USER:$CHEM_USER $VIRTUAL_ENV $VIRTUAL_ENV

ENV PATH=${VIRTUAL_ENV}/bin:$PATH
ENV CHEM_HOST='0.0.0.0'
ENV CHEM_PORT=45001
ENV CHEM_NUM_CORES=4

COPY --chown=$CHEM_USER:$CHEM_USER . .
RUN python -c "import rdkit.Chem.Draw, app"
EXPOSE 45001

ENTRYPOINT gunicorn app:app \
             --timeout 900 \
             -b $CHEM_HOST:$CHEM_PORT \
             --access-logfile - \
             --error-logfile - \
             --workers 2 \
             --threads=4 \
             --worker-class=gthread
