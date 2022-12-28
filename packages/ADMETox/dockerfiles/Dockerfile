FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml .
RUN conda env create -f environment.yml

SHELL ["conda","run","-n","myenv","/bin/bash","-c"]

COPY ADMETox.py .

#start server inside container
CMD ["conda", "run", "--no-capture-output","-n", "myenv", "python","ADMETox.py", "runserver","0.0.0.0:8000"]