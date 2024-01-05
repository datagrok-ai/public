# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04 AS builder

# Set noninteractive mode for apt-get
ENV DEBIAN_FRONTEND noninteractive

# Install required dependencies
RUN apt-get update && apt-get install -y \
    git \
    libopenbabel-dev \
    libopenbabel6 \
    pymol \
    python3-distutils \
    python3-lxml \
    python3-openbabel \
    python3-pymol \
    python3

# Install Flask
RUN apt-get install -y python3-pip
RUN pip3 install Flask

# Copy flask app
COPY plip.py src/plip.py

# Clone PLIP code from GitHub
RUN git clone https://github.com/pharmai/plip.git src/plip

# Set PYTHONPATH to include /src
ENV PYTHONPATH $PYTHONPATH:/src

# Expose port 8000
EXPOSE 8000

# Set entry point to run Flask app
ENTRYPOINT ["python3", "src/plip.py"]