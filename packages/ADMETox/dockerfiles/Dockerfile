FROM mambaorg/micromamba

SHELL ["/bin/bash","-c"]

USER root

RUN mkdir /app && chown mambauser:mambauser /app

USER mambauser

WORKDIR /app

RUN micromamba install -n base -c conda-forge awscli -y && \
    micromamba run -n base aws s3 cp s3://datagrok-data/models/admetox /app --no-sign-request --recursive && \
    micromamba remove -n base awscli -y && \
    micromamba clean -afy

COPY environment.yml /app/
COPY requirements.txt /app/

RUN set -ex ; \
    micromamba env create -f environment.yml ; \
    micromamba clean -afy

COPY . /app

SHELL ["micromamba","run","-n","myenv", "/bin/bash","-c"]

EXPOSE 8000

#ARG GUNICORN_MAX_REQUESTS=2
#ARG GUNICORN_MAX_REQUESTS_JITTER=3
ARG GUNICORN_WORKERS=4
ARG GUNICORN_THREADS=4
#ENV GUNICORN_MAX_REQUESTS="${GUNICORN_MAX_REQUESTS}"
#ENV GUNICORN_MAX_REQUESTS_JITTER="${GUNICORN_MAX_REQUESTS_JITTER}"
ENV GUNICORN_WORKERS="${GUNICORN_WORKERS}"
ENV GUNICORN_THREADS="${GUNICORN_THREADS}"

#start server inside container
ENTRYPOINT gunicorn ADMETox:application \
           --timeout 900 -b "0.0.0.0:8000" \
           --access-logfile - --error-logfile - \
           --workers=${GUNICORN_WORKERS} --threads=${GUNICORN_THREADS} --worker-class=gthread
#           --max-requests=${GUNICORN_MAX_REQUESTS} --max-requests-jitter=${GUNICORN_MAX_REQUESTS_JITTER}
