#=====================================================
#	FSM2oshd compilation script for Linux
#
# Simon Filhol, September 2023
# Compilation script based on:  https://github.com/RichardEssery/FSM/blob/master/compil.sh
# 
# This for the fsm2oshd version based on text files. Should be default. 
#
# requires:
# 	- gfortran
#=====================================================

FSM_NAME=FSM_OSHD

gfortran -o $FSM_NAME -O3 \
src/core/MODULES.F90 src/txt/MODULES.F90 \
src/core/LWRADTOPO.F90 src/core/SWRADTOPO.F90 src/core/CANOPY.F90 src/core/EBALFOR.F90 \
src/core/EBALSRF_SBG.F90 src/core/EBALSRF.F90 src/core/FSM2.F90 src/core/LUDCMP.F90 src/core/QSAT.F90 \
src/core/RADIATION.F90 src/core/SNOW.F90 src/core/SNOWCOVERFRACTION.F90 src/core/SOIL.F90 src/core/SFEXCH.F90 \
src/core/THERMAL.F90 src/core/TRIDIAG.F90 src/core/PHYSICS.F90 \
src/txt/OPEN_FILES.F90 src/txt/CUMULATE.F90 src/txt/DRIVE.F90 src/txt/DUMP.F90 src/txt/SETUP.F90 src/txt/OUTPUT.F90
rm *.mod
rm *.mod0

echo "---> FSM2oshd compiled as executable $FSM_NAME"
