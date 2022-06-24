FROM python:latest

RUN mkdir /app

COPY ./requirements.txt /app/requirements.txt
WORKDIR /app

RUN pip install -r requirements.txt


ADD src/* /app/src/
ADD data/* /app/data/

ENTRYPOINT ["python", "src/clustering_main.py", "--network", "data/DREAM_files/dream_3.txt"]