#!/bin/bash
#docker system prune -a
cd advdiff
docker build . -t actomtoolbox/adv-diff
cd ../carbonatesys
docker build . -t actomtoolbox/actom-co2
cd ../cseep
docker build . -t actomtoolbox/cseep
cd ../impacts
docker build . -t actomtoolbox/actom-impacts
cd ../optcover
docker build . -t actomtoolbox/opt-cover-app
cd ../output
docker build . -t actomtoolbox/actom-output
cd ../rerun-pyth
docker build . -t actomtoolbox/re-run
cd ../roc
docker build . -t actomtoolbox/actom-roc
cd ../setup
docker build . -t actomtoolbox/actom-run


docker push actomtoolbox/adv-diff
docker push actomtoolbox/actom-co2
docker push actomtoolbox/cseep
docker push actomtoolbox/actom-impacts
docker push actomtoolbox/opt-cover-app
docker push actomtoolbox/actom-output
docker push actomtoolbox/re-run
docker push actomtoolbox/actom-roc
docker push actomtoolbox/actom-run

