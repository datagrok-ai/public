FROM python:3.10-slim


RUN pip install pika
RUN pip install rdkit==2023.3.2
RUN apt-get update && apt-get install -y libxrender1 libxext6 libsm6 libglib2.0-0 libexpat1 \
    && rm -rf /var/lib/apt/lists/*
RUN pip install matplotlib\<3.8
RUN pip install numpy\<2
RUN pip install websockets
RUN pip install pandas
RUN pip install scikit-learn
COPY server.py server.py
EXPOSE 8080
CMD ["python3", "server.py"]
