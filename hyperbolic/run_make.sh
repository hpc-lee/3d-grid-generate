#!/bin/bash

#export GNUHOME=/data3/lihl/software/gcc-10.3.0-compile
#export NETCDF=/data3/lihl/software/disable-netcdf-4.4.1

export GNUHOME=/usr
export NETCDF=/data/apps/NetCDF/disable-netcdf-4.8.1

echo
echo "start to make ..."
make -j 
