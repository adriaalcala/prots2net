FROM python:3.8

RUN apt-get update && apt-get -y upgrade
RUN apt-get install libgirepository1.0-dev gcc libcairo2-dev pkg-config python3-dev gir1.2-gtk-3.0 -y
RUN apt-get install -y python3-gi
RUN pip install --no-cache-dir pipenv==2018.11.26

WORKDIR /app

COPY Pipfile Pipfile.lock ./
RUN pipenv install --system

RUN apt-get install -y gir1.2-webkit2-4.0
RUN pip3 install PyGObject

ENV PORT=8080
CMD python -m GUI.main

COPY . ./