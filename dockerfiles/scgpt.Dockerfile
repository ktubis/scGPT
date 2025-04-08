FROM xueerchen/scgpt:0.1.7
WORKDIR /scgpt
COPY . .
RUN pip install jupyter

