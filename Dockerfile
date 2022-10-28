# condaで環境構築
# deepTools
FROM continuumio/anaconda3:latest

MAINTAINER "edge2992"

# When use on ubuntu,
# UID & GID needs to be the same with your host environment
# to use git
ARG UNAME=app
ARG GROUPNAME=app
ARG UID=1000
ARG GID=1000
RUN groupadd -g $GID $GROUPNAME && \
  useradd -m -s /bin/bash -u $UID -g $GID $UNAME

RUN mkdir /opt/conda/envs/cq /opt/conda/pkgs && \
    chgrp $UNAME /opt/conda/pkgs && \
    chmod g+w /opt/conda/pkgs && \
    touch /opt/conda/pkgs/urls.txt && \
    chown $UNAME /opt/conda/pkgs/urls.txt

# COPY --chown=${UNAME}:${GROUPNAME} environment.yml labextensions.txt

RUN apt-get update && apt-get install -y vim wget sudo

# ENV PATH=/opt/conda/bin:${PATH}

RUN conda update -n base -c defaults conda

USER $UNAME
WORKDIR /home/$UNAME

# COPY protein.yml .
# RUN conda env create -n bio -f bio.yml

ARG CONDA_ENV=bio
RUN conda create --name $CONDA_ENV

RUN conda init bash
# RUN echo "conda activate ${CONDA_ENV}" >> ~/.bashrc


ENV CONDA_DEFAULT_ENV bio %% \
  PATH /home/${UNAME}/.conda/envs/${CONDA_ENV}/bin:$PATH

CMD ["/bin/bash"]
