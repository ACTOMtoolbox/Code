FROM python:3.7-slim


WORKDIR /app

COPY Requirements.txt /app/
COPY SetCoverForACTOM.py /app/
COPY Input /app/
COPY Output /app/

RUN pip install --no-cache-dir -r /app/Requirements.txt

CMD ["python","SetCoverForACTOM.py"]
