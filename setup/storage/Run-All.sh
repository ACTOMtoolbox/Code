echo '***************************************************'
echo '*                                                 *'
echo '*          Starting up the ACTOM Toolbox          *'
echo '*                   v1.0 2019-2022                *'
echo '*                                                 *'
echo '***************************************************'
echo
echo
#1
echo '***************************************************'
echo '*         The Advection-Diffusion Module          *'
echo '*          (Tracer Transport Simulator)           *'
echo '***************************************************'
echo
mkdir -p ../Advdiff
cd ../Advdiff
mv ../Setup.ini AdvDiff.ini
mkdir -p Figures
mkdir -p output
docker run -it $options \
          #insertttminputhere
          --mount type=bind,source="$(pwd)"/,target=/app/External-Indata \
          --mount type=bind,source="$(pwd)"/output,target=/app/Outdata \
          --mount type=bind,source="$(pwd)"/Figures,target=/app/Figures \
 actomtoolbox/adv-diff
echo
#2
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
#3
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
#4
echo '***************************************************'
echo '*              Deployment strategies              *'
echo '*                  Optimal Cover                  *'
echo '***************************************************'
echo
mkdir -p ../OptCover
cd ../OptCover
mkdir -p output
docker run -it $options \
          #insertoptinputhere
          --mount type=bind,source="$(pwd)"/output,target=/app/Output \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
 actomtoolbox/opt-cover-app
echo
#5
echo '***************************************************'
echo '*              The Carbonate System               *'
echo '*                                                 *'
echo '***************************************************'
echo
mkdir -p ../carbon
cd ../carbon
mkdir -p output
docker run -it $options \
          #insertcarbinputhere
          --mount type=bind,source="$(pwd)"/output,target=/external/output \
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
 actomtoolbox/actom-co2
echo
#6
echo '***************************************************'
echo '*                 Impact Analysis                 *'
echo '*                                                 *'
echo '***************************************************'
echo
mkdir -p ../impacts
cd ../impacts
mkdir -p output
docker run -it $options \
          #insertimpinputhere
          --mount type=bind,source="$(pwd)"/..,target=/external/settings \
          --mount type=bind,source="$(pwd)"/output,target=/external/output \
 actomtoolbox/actom-impacts
echo
#7
echo '***************************************************'
echo '*                 OUTPUT - Report                 *'
echo '*                                                 *'
echo '***************************************************'
echo
 cd ../
docker run -it $options -e systemOS=$systemOS --mount type=bind,source="$(pwd)",target=/srv/actom-output/input actomtoolbox/actom-output
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
echo "*                 'sh Re-Run.sh'                   *"
echo '*                                                 *'
echo '***************************************************'
echo
#8
