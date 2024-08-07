#!/bin/bash

#export MPIHOME=/data3/lihl/software/openmpi-gnu-4.1.2
#export NETCDF=/data3/lihl/software/disable-netcdf-4.4.1

export MPIHOME=/data/apps/openmpi/4.1.5-cuda-aware
export NETCDF=/data/apps/NetCDF/disable-netcdf-4.8.1

echo
echo "start to make ..."
make -j 
