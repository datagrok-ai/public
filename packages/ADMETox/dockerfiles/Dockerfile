# base image
FROM python:2.7-alpine

# setup environment variable  
ENV APP_HOME=/home/app/webapp  

# where your code lives  
WORKDIR $APP_HOME

# set work directory  
COPY ADMETox.py ./

COPY requirements.txt ./

# set environment variables  
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1  

# install dependencies 
RUN apk add build-base linux-headers --no-cache 
RUN python -m pip install --upgrade pip  

# run this command to install all dependencies  
RUN python -m pip install -r $APP_HOME/requirements.txt   

# start server  
ENTRYPOINT [ "python", "ADMETox.py", "runserver", "0.0.0.0:8000" ]

# port where the Django app runs
EXPOSE 8000