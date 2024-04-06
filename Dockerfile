# Start from the Jupyter base image
FROM jupyter/scipy-notebook:latest

USER root

# Install Java
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk curl unzip wget nano && \
    apt-get clean;
    
RUN pip install rdkit

ARG BIOVER="07d5e563edb0"

RUN wget https://bitbucket.org/wishartlab/biotransformer3.0jar/get/${BIOVER}.zip && \
    unzip ${BIOVER}.zip && \
    rm ${BIOVER}.zip 

RUN mv wishartlab-biotransformer3.0jar-${BIOVER} mybio && \
    cd mybio

# Switch back to the jovyan user
USER jovyan
