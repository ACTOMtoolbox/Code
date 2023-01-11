echo '***************************************************'
echo '*                                                 *'
echo '*          Starting up the ACTOM Toolbox          *'
echo '*                   v1.0 2019-2022                *'
echo '*                                                 *'
echo '***************************************************'
echo
echo
echo '***************************************************'
echo '*         The Advection-Diffusion Module          *'
echo '*          (Tracer Transport Simulator)           *'
echo '***************************************************'
echo
mkdir -p Advdiff
cd Advdiff
cp ../../AdvDiff.ini AdvDiff.ini
mkdir -p Figures
mkdir -p output
docker run -it $options \
          --mount type=bind,source="$(pwd)",target=/external/settings \
          --mount type=bind,source="$(pwd)"/output,target=/app/Outdata \
          --mount type=bind,source="$(pwd)"/Figures,target=/app/Figures \
 actomtoolbox/adv-diff
echo
echo '***************************************************'
echo '*         Anomaly Criteria Identification         *'
echo '*                      Cseep                      *'
echo '***************************************************'
echo
mkdir -p ../Cseep
cd ../Cseep
mkdir -p output
docker run -it $options \
          #insertcseepinputhere
          --mount type=bind,source="$(pwd)"/output,target=/srv/actom-app/output \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
 actomtoolbox/cseep
echo
echo '***************************************************'
echo '*         Anomaly Criteria Identification         *'
echo '*                  Rate of Change                 *'
echo '***************************************************'
echo
mkdir -p ../ROC
cd ../ROC
mkdir -p output
docker run -it $options \
          #insertrocinputhere
          --mount type=bind,source="$(pwd)"/output,target=/external/output \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
 actomtoolbox/actom-roc
echo
echo '***************************************************'
echo '*              Deployment strategies              *'
echo '*                  Optimal Cover                  *'
echo '***************************************************'
echo
mkdir -p ../OptCover
cd ../OptCover
mkdir -p output
docker run -it $options \
          --mount type=bind,source="$(pwd)"/../Advdiff/output,target=/app/Input  \
          --mount type=bind,source="$(pwd)"/output,target=/app/Output \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
 actomtoolbox/opt-cover-app
echo
echo '***************************************************'
echo '*              The Carbonate System               *'
echo '*                                                 *'
echo '***************************************************'
echo
mkdir -p ../carbon
cd ../carbon
mkdir -p output
docker run -it $options \
          --mount type=bind,source="$(pwd)"/../Advdiff/output,target=/external/input \
          --mount type=bind,source="$(pwd)"/output,target=/external/output \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
 actomtoolbox/actom-co2
echo
echo '***************************************************'
echo '*                 Impact Analysis                 *'
echo '*                                                 *'
echo '***************************************************'
echo
mkdir -p ../impacts
cd ../impacts
mkdir -p output
docker run -it $options \
          --mount type=bind,source="$(pwd)"/../carbon/output,target=/external/input \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
          --mount type=bind,source="$(pwd)"/output,target=/external/output \
 actomtoolbox/actom-impacts
echo
echo '***************************************************'
echo '*                 OUTPUT - Report                 *'
echo '*                                                 *'
echo '***************************************************'
echo
 cd ../
 docker run -it $options --mount type=bind,source="$(pwd)",target=/srv/actom-output/input actomtoolbox/actom-output
echo
echo '***************************************************'
echo '*     The Technical Summary can be found at:      *'
echo '*                                                 *'
echo '* file://'$(pwd)'/Technical-Summary.html'
echo '*                                                 *'
echo '*         with the full report found at:          *'
echo '*                                                 *'
echo '* file://'$(pwd)'/Report.html'
echo '*                                                 *'
echo '*          press Ctrl and click on link           *'
echo '*     or Copy url into your favourite browser     *'
echo '*                                                 *'
echo '*                                                 *'
echo '*      To re-run the toolbox at a different       *'
echo '*  leakage rate, or with different pH thresholds  *'
echo '*             run the scripts below:              *'
echo "* 'cd $(pwd)' "
echo "*                 'sh Re-Run.sh'                  "
echo '*                                                 *'
echo '***************************************************'
echo
