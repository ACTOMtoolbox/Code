FROM python:3.8-slim
# surpress terminal error messages
ENV DEBIAN_FRONTEND noninteractive
# update the os
RUN apt-get clean \
    && apt-get -y update
# grab python packages that we need
RUN apt-get -y install python3-dev
RUN apt-get -y install python3-pip
RUN apt-get -y install python3-tk
RUN pip3 install Jinja2
RUN pip3 install pillow
RUN pip3 install tkhtmlview
#set default directory
WORKDIR /srv/actom-run
RUN chmod 777 /srv/actom-run
# copy codebase over
COPY . /srv/actom-run/
# set perms
RUN chmod +x start.sh
CMD ["./start.sh"]
