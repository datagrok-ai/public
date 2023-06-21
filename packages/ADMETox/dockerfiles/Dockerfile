FROM continuumio/miniconda3

WORKDIR /app

SHELL ["/bin/bash","-c"]

RUN pip install --upgrade awscli && \
    aws s3 cp s3://datagrok-data/models/admetox /app --no-sign-request --recursive && \
    pip uninstall -y awscli

COPY environment.yml /app/
COPY requirements.txt /app/

RUN set -ex ; \
    conda install -n base conda-libmamba-solver -y ; \
    conda config --set solver libmamba ; \
    conda env create -f environment.yml ; \
    conda clean -afy

COPY . /app

SHELL ["conda","run","-n","myenv","/bin/bash","-c"]

EXPOSE 8000

#start server inside container
ENTRYPOINT ["conda","run","--no-capture-output","-n","myenv","gunicorn", "ADMETox:application","--timeout","900","-b", "0.0.0.0:8000","--access-logfile", "-","--error-logfile", "-","--workers=2","--threads=4","--worker-class=gthread","--max-requests=100","--max-requests-jitter=15"]
