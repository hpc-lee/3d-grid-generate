#!/bin/bash

#set -x
set -e

date

#-- program related dir
EXEC_GRID=`pwd`/../main
echo "EXEC_GRID=${EXEC_GRID}"

#-- input dir
INPUTDIR1=/data/lihl/code/3d-grid-generate/elliptic/project1/output
INPUTDIR2=/data/lihl/code/3d-grid-generate/elliptic/project2/output

#-- output and conf
PROJDIR=`pwd`/../project
PAR_FILE=${PROJDIR}/test.json
OUTPUT_DIR=${PROJDIR}/output

rm -rf ${PROJDIR}

#-- create dir
mkdir -p ${PROJDIR}
mkdir -p ${OUTPUT_DIR}

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "input_grids_info" : [
    {
      "grid_import_dir" : "${INPUTDIR1}",
      "number_of_grid_points" : [250,250,125],
      "number_of_mpiprocs_in" : [2,2,2]
    },
    {
      "grid_import_dir" : "${INPUTDIR2}",
      "number_of_grid_points" : [250,250,200],
      "number_of_mpiprocs_in" : [3,2,2]
    }
  ],
    
  "merge_direction" : "z",

  "number_of_mpiprocs_out" : [3,3,2],

  "pml_layers" : {
         "number_of_pml_x1" : 10,
         "number_of_pml_x2" : 10,
         "number_of_pml_y1" : 10,
         "number_of_pml_y2" : 10,
         "number_of_pml_z1" : 10,
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

  "flag_sample" : 1,
  "sample_factor_xi" : 2,
  "sample_factor_et" : 2,
  "sample_factor_zt" : 3,

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
