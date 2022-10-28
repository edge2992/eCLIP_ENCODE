# condaで環境構築
# deepTools
FROM continuumio/anaconda3:latest

MAINTAINER "edge2992"

RUN apt-get update && apt-get install -y vim wget

WORKDIR /opt

ENV PATH=/opt/conda/bin:${PATH}

RUN conda update -n base -c defaults conda

# COPY protein.yml .
# RUN conda env create -n bio -f bio.yml

RUN conda create --name bio

RUN  conda init && \
  echo "conda activate bio" >> ~/.bashrc


ENV CONDA_DEFAULT_ENV bio %% \
  PATH /opt/conda/envs/bio/bin:$PATH

CMD ["/bin/bash"]
