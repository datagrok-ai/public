FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml .
COPY requirements.txt .
RUN set -ex ; \
    conda install -n base conda-libmamba-solver -y ; \
    conda config --set solver libmamba ; \
    conda env create -f environment.yml ; \
    conda clean -afy

SHELL ["conda","run","-n","myenv","/bin/bash","-c"]

COPY ADMETox.py .

#start server inside container
CMD ["conda", "run", "--no-capture-output","-n", "myenv", "python","ADMETox.py", "runserver","0.0.0.0:8000"]