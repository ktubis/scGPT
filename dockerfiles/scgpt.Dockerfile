FROM nvcr.io/nvidia/pytorch:23.05-py3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
    git \
    r-base \
    tzdata \
    build-essential

RUN pip install --upgrade pip
RUN pip install scgpt "flash-attn<1.0.5" "orbax<0.1.8" "numpy<2"
RUN pip install torch==2.2.0 torchtext==0.17.0 --upgrade
RUN pip install wandb
RUN pip install git+https://github.com/thunlp/OpenDelta.git
COPY tmp_delta_basemodel.py /tmp/tmp_delta_basemodel.py
RUN cp /tmp/tmp_delta_basemodel.py /usr/local/lib/python3.10/dist-packages/opendelta/basemodel.py


WORKDIR /scgpt


#FROM xueerchen/scgpt:0.1.7

#RUN pip install jupyter

