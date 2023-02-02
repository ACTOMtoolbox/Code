@ robocopy  ..\Advdiff Advdiff /E /DCOPY:DAT >NUL
@ robocopy  ..\ROC ROC /E /DCOPY:DAT >NUL
@ robocopy  ..\Cseep Cseep /E /DCOPY:DAT >NUL>NUL
@ echo.
@ echo ***************************************************
@ echo *                                                 *
@ echo *          Starting up the ACTOM Toolbox          *
@ echo *                      Re-Run                     *
@ echo *                   v1.0 2019-2022                *
@ echo *                                                 *
@ echo *         The Advection-Diffusion Module          *
@ echo *           (Tracer Transport Simulator)          *
@ echo *         Anomaly Criteria Identification         *
@ echo *           (Cseep and Rate of Change)            *
@ echo *            All give the same results            *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ echo.
@ echo ***************************************************
@ echo *              Deployment strategies              *
@ echo *                  Optimal Cover                  *
@ echo ***************************************************
@ echo.
@ mkdir OptCover
@ cd OptCover
@ mkdir output
docker run -i $options --mount type=bind,source="%cd%"/../Advdiff/output,target=/app/Input --mount type=bind,source=%cd%/output,target=/app/Output --mount type=bind,source=%cd%/..,target=/external/settings actomtoolbox/opt-cover-app
@ echo.
@ echo ***************************************************
@ echo *              The Carbonate System               *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ mkdir ..\carbon
@ cd ../carbon
@ mkdir output
docker run -i $options --mount type=bind,source="%cd%"/../Advdiff/output,target=/external/input --mount type=bind,source=%cd%/output,target=/external/output --mount type=bind,source=%cd%/..,target=/external/settings actomtoolbox/actom-co2
@ echo.
@ echo ***************************************************
@ echo *                 Impact Analysis                 *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ mkdir ..\impacts
@ cd ../impacts
@ mkdir output
docker run -i $options --mount type=bind,source="%cd%"/../carbon/output,target=/external/input --mount type=bind,source=%cd%/..,target=/external/settings --mount type=bind,source=%cd%/output,target=/external/output actomtoolbox/actom-impacts
@ echo.
@ echo ***************************************************
@ echo *                 OUTPUT - Report                 *
@ echo *                                                 *
@ echo ***************************************************
@ echo.
@ cd ../
docker run -i $options --mount type=bind,source=%cd%,target=/srv/actom-output/input actomtoolbox/actom-output
@ DEL Re-Run.sh
@ echo.
@ echo ***************************************************
@ echo *     The Technical Summary can be found at:      *
@ echo *                                                 *
@ echo * file://%cd%/Technical-Summary.html'
@ echo *                                                 *
@ echo *         with the full report found at:          *
@ echo *                                                 *
@ echo * file://%cd%/Report.html'
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
