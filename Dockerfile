FROM python:3.10.4-slim AS compile-image
ENV TZ America/Sao_Paulo
RUN apt-get update
RUN apt-get install -y build-essential gcc libbz2-dev libcurl4-openssl-dev zlib1g-dev liblzma-dev

RUN python -m venv /opt/venv
# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

COPY requirements.txt /requirements.txt
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r /requirements.txt

FROM python:3.10.4-slim

COPY --from=compile-image /opt/venv /opt/venv

COPY make_diploid.py /usr/bin/make_diploid.py
RUN chmod +x /usr/bin/make_diploid.py

ENV PATH="/opt/venv/bin:$PATH"
