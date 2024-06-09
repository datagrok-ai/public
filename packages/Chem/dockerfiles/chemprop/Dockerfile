FROM mambaorg/micromamba

SHELL ["/bin/bash", "-c"]

USER root

RUN mkdir /app && chown mambauser:mambauser /app

USER mambauser

WORKDIR /app

COPY --chown=mambauser:mambauser environment.yml /app/


RUN micromamba install -n base -c conda-forge awscli -y && \
    micromamba remove -n base awscli -y && \
    micromamba clean -afy

RUN set -ex ; \
    micromamba env create -f environment.yml ; \
    micromamba clean -afy

COPY --chown=mambauser:mambauser . /app

SHELL ["micromamba", "run", "-n", "myenv", "/bin/bash", "-c"]

EXPOSE 5000

ENTRYPOINT ["micromamba", "run", "-n", "myenv", "python", "grok_compute.py"]
