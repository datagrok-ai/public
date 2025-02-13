FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    gcc \
    git \
    python3-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN python -m venv /env

RUN /env/bin/python -m pip install --upgrade pip \
    && git clone https://github.com/jwohlwend/boltz.git \
    && cd boltz \
    && /env/bin/pip install -e .

RUN /env/bin/pip install flask
RUN /env/bin/pip install Flask-Cors

COPY . /app

WORKDIR /app

ENV PATH="/env/bin:$PATH"

EXPOSE 8001

ENTRYPOINT ["python", "boltz.py", "--host", "0.0.0.0", "--port", "8001"]