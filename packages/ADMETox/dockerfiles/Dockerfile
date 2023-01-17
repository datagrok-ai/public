FROM continuumio/miniconda3

WORKDIR /app

COPY . /app

RUN set -ex ; \
    conda install -n base conda-libmamba-solver -y ; \
    conda config --set solver libmamba ; \
    conda env create -f environment.yml ; \
    conda clean -afy

SHELL ["conda","run","-n","myenv","/bin/bash","-c"]

EXPOSE 8000

#start server inside container
ENTRYPOINT ["conda", "run", "--no-capture-output","-n", "myenv", "python","ADMETox.py", "runserver","0.0.0.0:8000"]
