FROM mambaorg/micromamba

SHELL ["/bin/bash","-c"]

USER root

RUN mkdir /app && chown mambauser:mambauser /app

USER mambauser

WORKDIR /app

COPY environment.yml /app/
COPY requirements.txt /app/

RUN micromamba install -n base -c conda-forge awscli -y && \
    micromamba run -n base aws s3 cp s3://datagrok-data/models/admetox /app --no-sign-request --recursive --region us-east-2 && \
    micromamba remove -n base awscli -y && \
    micromamba clean -afy

RUN set -ex ; \
    micromamba env create -f environment.yml ; \
    micromamba clean -afy

COPY . /app

SHELL ["micromamba","run","-n","myenv", "/bin/bash","-c"]

EXPOSE 8000

ENTRYPOINT ["micromamba", "run", "-n", "myenv", "gunicorn", "-w", "1", "-b", "0.0.0.0:8000", "--timeout", "3000", "ADMETox:app"]