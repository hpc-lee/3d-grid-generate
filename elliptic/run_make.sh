#!/bin/bash

export MPIHOME=/data/apps/openmpi/4.1.5-cuda-aware
export NETCDF=/data/apps/NetCDF/disable-netcdf-4.8.1

echo
echo "start to make ..."
make -j 
