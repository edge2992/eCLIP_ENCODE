FROM continuumio/miniconda3:4.12.0
MAINTAINER "edge2992"

# Even if you change UID or GID later,
# only the ownership under /home/cq will change, so define it here.
ARG UID=1000
ARG GID=1000
RUN groupadd -g $GID cq && \
  adduser --home /home/cq --disabled-password \
  --gecos "Default user" --ingroup  cq --uid $UID cq

RUN apt-get update && apt-get -y install gosu sudo

RUN conda update -n base -c defaults conda

RUN mkdir /opt/conda/envs/cq && \
    chgrp -R cq /opt/conda/pkgs && \
    chmod -R g+w /opt/conda/pkgs && \
    touch /opt/conda/pkgs/urls.txt && \
    chown cq /opt/conda/envs/cq /opt/conda/pkgs/urls.txt

COPY --chown=cq:cq environment.yml labextensions.txt /opt/conda/envs/cq/

WORKDIR /home/cq
USER cq

RUN conda env create -f /opt/conda/envs/cq/environment.yml -n cq

# instead of running "RUN conda init bash"
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate cq" >> ~/.bashrc

ENV CONDA_DEFAULT_ENV cq %% \
  PATH /opt/conda/envs/cq/bin:$PATH

COPY ./entrypoint.sh /tmp

ENTRYPOINT ["/tmp/entrypoint.sh"]

CMD ["/bin/bash"]
