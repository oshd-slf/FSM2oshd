::-------------------------------------------------------------------------------------------
:: Compile using the same setup as in jim_operational...
::-------------------------------------------------------------------------------------------
echo off
set optim=-O3

echo on
set mods= core/MODULES.F90 oshd/MODULES.F90
set routines= core/LWRADTOPO.F90 core/SWRADTOPO.F90 core/CANOPY.F90 core/EBALFOR.F90 ^
core/EBALSRF_SBG.F90 core/EBALSRF.F90 core/FSM2.F90 core/LUDCMP.F90 core/QSAT.F90 ^
core/RADIATION.F90 core/SNOW.F90 core/SNOWCOVERFRACTION.F90 core/SOIL.F90 core/SFEXCH.F90 ^
core/THERMAL.F90 core/TRIDIAG.F90 core/PHYSICS.F90 ^
oshd/OPEN_FILES.F90 oshd/CUMULATE.F90 oshd/DRIVE.F90 oshd/DUMP.F90 oshd/SETUP.F90 oshd/OUTPUT.F90
if /I "%1" == "-O3" set optim= %1
if /I "%1" == "-O2" set optim= %1
if /I "%1" == "-O1" set optim= %1
if /I "%1" == "-O0" set optim= %1
set profil= %2
cd src
gfortran %mods% %routines% %optim% %profil% -cpp -ffpe-trap=overflow -o FSM_OSHD
del *.mod
move FSM_OSHD.exe ..\FSM_OSHD.exe
cd ..