FROM python:3.7-slim

# surpress terminal error messages
ENV DEBIAN_FRONTEND noninteractive

# update the os
RUN apt-get clean \
    && apt-get -y update

# grab python packages that we need
RUN apt-get -y install python3-dev

# COPY python/PIP requirements over in a text file
COPY requirements.txt /srv/actom-app/requirements.txt
# install python libraries
RUN pip install -r /srv/actom-app/requirements.txt --src /usr/local/src

#set default directory
WORKDIR /srv/actom-app

# copy codebase over
COPY . /srv/actom-app/

# set perms
RUN chmod +x start.sh
CMD ["./start.sh"]
