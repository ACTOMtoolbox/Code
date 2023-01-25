# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3.8.12

# Keeps Python from generating .pyc files in the container
#ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
#ENV PYTHONUNBUFFERED=1

WORKDIR /app
COPY requirements.txt .

ENV DEBIAN_FRONTEND noninteractive
RUN python -m pip install --upgrade pip \
&& python -m pip install --upgrade setuptools \
&& apt-get update -y \
&& apt-get install build-essential \
&& apt install libgl1-mesa-glx -y \
&& apt-get install 'ffmpeg' \
    'libsm6' \
    'libxext6'  -y \
&& apt-get install libglu1 -y \
&& cat requirements.txt | xargs -n 1 -L 1 pip install


# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
#RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
#USER appuser
COPY . /app
# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
CMD ["python", "main.py"]

