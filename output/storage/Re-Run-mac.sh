#!/bin/bash
TIMESTAMP=$(date +%Y%m%d%H%M%S)
title=re-run-$TIMESTAMP

mkdir -p $title
cp data.ini $title/data.ini
cp -r logs/ $title/logs
cd $title

IP=$(/usr/sbin/ipconfig getifaddr en0)
/opt/X11/bin/xhost + "$IP"

docker run -i -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY="${IP}:0" -e systemOS=mac --mount type=bind,source="$(pwd)",target=/srv/actom-rerun/input actomtoolbox/re-run

sh $(pwd)/Re-Run-Toolbox.sh
mv $(pwd)/Re-Run-Toolbox.sh logs/Re-Run-Toolbox.sh
