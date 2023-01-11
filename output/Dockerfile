FROM python:3.7-slim

# surpress terminal error messages
ENV DEBIAN_FRONTEND noninteractive

# update the os
RUN apt-get clean \
    && apt-get -y update

# grab python packages that we need
RUN apt-get -y install python3-dev
RUN apt-get -y install python3-pip
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install Jinja2
RUN pip3 install matplotlib
RUN pip3 install netCDF4
RUN pip3 install termcolor

#set default directory
WORKDIR /srv/actom-output

RUN chmod 777 /srv/actom-output

# copy codebase over
COPY . /srv/actom-output/

# set perms
RUN chmod +x start.sh
CMD ["./start.sh"]
