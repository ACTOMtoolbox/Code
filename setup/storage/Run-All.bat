@ echo ***************************************************
@ echo *                                                 *
@ echo *          Starting up the ACTOM Toolbox          *
@ echo *                   v1.0 2019-2022                *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ echo.
#1
@ echo ***************************************************
@ echo *         The Advection-Diffusion Module          *
@ echo *          (Tracer Transport Simulator)           *
@ echo ***************************************************
@ echo.
@ mkdir ../Advdiff
@ move Setup.ini %cd%/Advdiff > nul
@ cd ../Advdiff
@ RENAME "Setup.ini" "AdvDiff.ini"
@ mkdir Figures
@ mkdir output
docker run -i $options #insertttminputhere --mount type=bind,source=%cd%/,target=/app/External-Indata --mount type=bind,source=%cd%/output,target=/app/Outdata --mount type=bind,source=%cd%/Figures,target=/app/Figures actomtoolbox/adv-diff
@ echo.
#2
@ echo ***************************************************
@ echo *         Anomaly Criteria Identification         *
@ echo *                      Cseep                      *
@ echo ***************************************************
@ echo.
@ mkdir ..\Cseep
@ cd ../Cseep
@ mkdir output
docker run -it $options #insertcseepinputhere --mount type=bind,source=%cd%/output,target=/srv/actom-app/output --mount type=bind,source=%cd%/..,target=/external/settings actomtoolbox/cseep
@ echo.
#3
@ echo ***************************************************
@ echo *         Anomaly Criteria Identification         *
@ echo *                  Rate of Change                 *
@ echo ***************************************************
@ echo.
@ mkdir ..\ROC
@ cd ../ROC
@ mkdir output
docker run -it $options #insertrocinputhere --mount type=bind,source=%cd%/output,target=/external/output --mount type=bind,source=%cd%/..,target=/external/settings actomtoolbox/actom-roc
@ echo.
#4
@ echo ***************************************************
@ echo *              Deployment strategies              *
@ echo *                  Optimal Cover                  *
@ echo ***************************************************
@ echo.
@ mkdir ..\OptCover
@ cd ../OptCover
@ mkdir output
docker run -it $options #insertoptinputhere --mount type=bind,source=%cd%/output,target=/app/Output --mount type=bind,source=%cd%/..,target=/external/settings actomtoolbox/opt-cover-app
@ echo.
#5
@ echo ***************************************************
@ echo *              The Carbonate System               *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ mkdir ..\carbon
@ cd ../carbon
@ mkdir output
docker run -it $options #insertcarbinputhere --mount type=bind,source=%cd%/output,target=/external/output --mount type=bind,source=%cd%/..,target=/external/settings actomtoolbox/actom-co2
@ echo.
#6
@ echo ***************************************************
@ echo *                 Impact Analysis                 *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ mkdir ..\impacts
@ cd ../impacts
@ mkdir output
docker run -it $options #insertimpinputhere --mount type=bind,source=%cd%/..,target=/external/settings --mount type=bind,source=%cd%/output,target=/external/output actomtoolbox/actom-impacts
@ echo.
#7
@ echo ***************************************************
@ echo *                 OUTPUT - Report                 *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ cd ../
docker run -it $options -e systemOS=$systemOS --mount type=bind,source=%cd%,target=/srv/actom-output/input actomtoolbox/actom-output
@ echo.
@ echo ***************************************************
@ echo *     The Technical Summary can be found at:      *
@ echo *                                                 *
@ echo * file://%cd%/Technical-Summary.html
@ echo *                                                 *
@ echo *         with the full report found at:          *
@ echo *                                                 *
@ echo * file://%cd%/Report.html
@ echo *                                                 *
@ echo *          press Ctrl and click on link           *
@ echo *     or Copy url into your favourite browser     *
@ echo *                                                 *
@ echo *                                                 *
@ echo *      To re-run the toolbox at a different       *
@ echo *  leakage rate, or with different pH thresholds  *
@ echo *             run the scripts below:              *
@ echo * 'cd %cd%'
@ echo *                 '.\Re-Run.bat'                  *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
#8
@ mkdir logs
@ copy Run-All.bat logs >NUL
@ del "%~f0" & exit
