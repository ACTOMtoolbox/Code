#!/bin/bash
TIMESTAMP=$(date +%Y%m%d%H%M%S)
title=re-run-$TIMESTAMP

mkdir -p $title
mkdir -p $title/logs
cp data.ini $title/data.ini
cp logs/Run-All.sh $title/logs/Run-All.sh
cd $title

docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$DISPLAY --mount type=bind,source="$(pwd)",target=/srv/actom-rerun/input actomtoolbox/re-run

sh $(pwd)/Re-Run-Toolbox.sh

rm -r Re-Run-Toolbox.sh
