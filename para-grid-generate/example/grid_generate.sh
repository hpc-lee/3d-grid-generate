#!/bin/bash

#set -x
set -e

date

#-- system related dir
MPIDIR=/data3/lihl/software/openmpi-gnu-4.1.2

#-- program related dir
EXEC_GRID=`pwd`/../main_grid_3d
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

#-- total x mpi procs
NPROCS_X=3
#-- total y mpi procs
NPROCS_Y=2
#-- total z mpi procs
NPROCS_Z=3
#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "number_of_grid_points_x" : 200,
  "number_of_grid_points_y" : 200,
  "number_of_grid_points_z" : 100,

  "number_of_mpiprocs_x" : $NPROCS_X,
  "number_of_mpiprocs_y" : $NPROCS_Y,
  "number_of_mpiprocs_z" : $NPROCS_Z,

  "#pml_layers" : {
         "number_of_pml_x1" : 20,
         "number_of_pml_x2" : 20,
         "number_of_pml_y1" : 20,
         "number_of_pml_y2" : 20,
         "number_of_pml_z1" : 20,
         "number_of_pml_z2" : 20
  },

  "check_orth" : 1,
  "check_jac" : 1,
  "check_step_xi" : 1,
  "check_step_et" : 1,
  "check_step_zt" : 1,
  "check_smooth_xi" : 1,
  "check_smooth_et" : 1,
  "check_smooth_zt" : 1,

  "geometry_input_file" : "${INPUTDIR}/data_file_3d.txt",
  "grid_export_dir" : "${OUTPUT_DIR}",

  "grid_method" : {
      "elli_diri" : {
          "coef" : -10,
          "iter_err" : 1E-2,
          "max_iter" : 5E3
      },
      "#elli_higen" : {
          "coef" : -20,
          "iter_err" : 1E-2,
          "max_iter" : 5E3,
          "distance" : [100,100,100,100,100,100]
      }
  }
}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#
#-- get np
NUMPROCS_X=`grep number_of_mpiprocs_x ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS_Y=`grep number_of_mpiprocs_y ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS_Z=`grep number_of_mpiprocs_z ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS=$(( NUMPROCS_X*NUMPROCS_Y*NUMPROCS_Z ))
echo $NUMPROCS_X $NUMPROCS_Y $NUMPROCS_Z $NUMPROCS

#-- gen run script
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