@ for /f "delims=[] tokens=2" %%a in ('ping -4 -n 1 %ComputerName% ^| findstr [') do @ set NetworkIP=%%a
@ set Disp="%NetworkIP%:0.0"

@ set TIMESTAMP=%DATE:~6,4%_%DATE:~3,2%_%DATE:~0,2%__%TIME:~0,2%_%TIME:~3,2%_%TIME:~6,2%
@ set TIMESTAMP=%TIMESTAMP: =0%
@ set TIMESTAMP=%TIMESTAMP:_=%
@ set title=re-run-%TIMESTAMP%

@ mkdir %title%
@ copy data.ini %title%\data.ini>NUL
@ robocopy logs\ %title%\logs>NUL
@ cd %title%

docker run -i -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=%Disp% --mount type=bind,source=%cd%,target=/srv/actom-rerun/input actomtoolbox/re-run

@ .\Re-Run-Toolbox.bat
@ move Run-All.bat %cd%/logs >NUL
