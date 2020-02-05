FROM continuumio/miniconda3

COPY . /app/

RUN conda env update --name base --file /app/environment.yml

ENTRYPOINT [ "/app/trimviz.py" ]