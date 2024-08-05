#!/bin/bash

#set -x
set -e

date

#-- program related dir
EXEC_GRID=`pwd`/../grid_post_proc_3d
echo "EXEC_GRID=${EXEC_GRID}"

#-- input dir
#INPUTDIR=/data3/lihl/code-lihl/3d-grid-generate/ellipitic/project/output
INPUTDIR=/data3/lihl/code-lihl/3d-grid-generate/parabolic/project/output

#-- output and conf
PROJDIR=`pwd`/../project
PAR_FILE=${PROJDIR}/test.json
OUTPUT_DIR=${PROJDIR}/output

rm -rf ${PROJDIR}

#-- create dir
mkdir -p ${PROJDIR}
mkdir -p ${OUTPUT_DIR}

# grid generate procs
#-- total x mpi procs
NPROCS_X_IN=1
#-- total y mpi procs
NPROCS_Y_IN=10
#-- total z mpi procs
NPROCS_Z_IN=1

# after post procs
#-- total x mpi procs
NPROCS_X_OUT=2
#-- total y mpi procs
NPROCS_Y_OUT=2
#-- total z mpi procs
NPROCS_Z_OUT=1

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "number_of_grid_points_x" : 400,
  "number_of_grid_points_y" : 300,
  "number_of_grid_points_z" : 200,

  "number_of_mpiprocs_x_in" : $NPROCS_X_IN,
  "number_of_mpiprocs_y_in" : $NPROCS_Y_IN,
  "number_of_mpiprocs_z_in" : $NPROCS_Z_IN,

  "number_of_mpiprocs_x_out" : $NPROCS_X_OUT,
  "number_of_mpiprocs_y_out" : $NPROCS_Y_OUT,
  "number_of_mpiprocs_z_out" : $NPROCS_Z_OUT,

  "pml_layers" : {
         "number_of_pml_x1" : 0,
         "number_of_pml_x2" : 0,
         "number_of_pml_y1" : 0,
         "number_of_pml_y2" : 0,
         "number_of_pml_z1" : 0,
         "number_of_pml_z2" : 0
  },

  "check_orth" : 1,
  "check_jac" : 1,
  "check_step_xi" : 1,
  "check_step_et" : 1,
  "check_step_zt" : 1,
  "check_smooth_xi" : 1,
  "check_smooth_et" : 1,
  "check_smooth_zt" : 1,


  "flag_strech_xi" : 0,
  "strech_xi_coef" : 0.0001,
  "flag_strech_et" : 0,
  "strech_et_coef" : 0.0001,
  "flag_strech_zt" : 0,
  "strech_zt_coef" : 0.0001,

  "flag_sample" : 0,
  "sample_factor_xi" : 1,
  "sample_factor_et" : 2,
  "sample_factor_zt" : 3,

  "grid_import_dir" : "${INPUTDIR}",
  "grid_export_dir" : "${OUTPUT_DIR}"

}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#
#-- gen run script
cat << ieof > ${PROJDIR}/grid_generate.sh
#!/bin/bash

set -e

printf "\nStart grid post process ...\n";
time $EXEC_GRID $PAR_FILE 100 2>&1 |tee log
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
