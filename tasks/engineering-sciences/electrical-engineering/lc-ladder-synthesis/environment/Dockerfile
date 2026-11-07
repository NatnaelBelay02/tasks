FROM python:3.11-slim

RUN pip install --no-cache-dir \
    numpy==2.2.6 \
    scipy==1.15.3 \
    sympy==1.13.3 \
    mpmath==1.3.0

WORKDIR /root

COPY data /root/data

RUN mkdir -p /root/results
