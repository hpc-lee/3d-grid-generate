#!/bin/bash

#set -x
set -e

date

#-- system related dir
#MPIDIR=/data3/lihl/software/openmpi-gnu-4.1.2
MPIDIR=/data/apps/openmpi/4.1.5-cuda-aware

#-- program related dir
EXEC_GRID=`pwd`/../main
echo "EXEC_GRID=${EXEC_GRID}"

#-- input dir
INPUTDIR=`pwd`

#-- output and conf
PROJDIR=`pwd`/../project
PAR_FILE=${PROJDIR}/test.json
OUTPUT_DIR=${PROJDIR}/output

rm -rf ${PROJDIR}

#-- create dir
mkdir -p ${PROJDIR}
mkdir -p ${OUTPUT_DIR}

#-- total mpi procs
NUMPROCS=1
#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "number_of_grid_points_x" : 1000,
  "number_of_grid_points_y" : 500,
  "number_of_grid_points_z" : 500,

  "number_of_mpiprocs" : $NUMPROCS,

  "check_orth" : 1,
  "check_jac"  : 0,
  "check_step_xi" : 0,
  "check_step_et" : 0,
  "check_step_zt" : 0,
  "check_smooth_xi" : 1,
  "check_smooth_et" : 1,
  "check_smooth_zt" : 1,

  "geometry_input_file" : "${INPUTDIR}/data_file_3d.txt",
  "step_input_file" : "${INPUTDIR}/step_file_3d.txt",
  "grid_export_dir" : "${OUTPUT_DIR}",

  "parabolic" : {
      "coef" : 60,
      "t2b" : 1,
      "direction" : "z"
  }
}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#
#-- gen run script
echo $NUMPROCS

cat << ieof > ${PROJDIR}/grid_generate.sh
#!/bin/bash

set -e

printf "\nUse $NUMPROCS CPUs on following nodes:\n"

printf "\nStart grid generate ...\n";
time $MPIDIR/bin/mpiexec -np $NUMPROCS $EXEC_GRID $PAR_FILE 100 2 2>&1 |tee log
if [ $? -ne 0 ]; then
    printf "\ngrid generate fail! stop!\n"
    exit 1
fi

ieof

#-------------------------------------------------------------------------------
#-- start run
#-------------------------------------------------------------------------------

chmod 755 ${PROJDIR}/grid_generate.sh
${PROJDIR}/grid_generate.sh

date

# vim:ts=4:sw=4:nu:et:ai:
