#!/bin/bash

#set -x
set -e

date

#-- program related dir
EXEC_GRID=`pwd`/../main
echo "EXEC_GRID=${EXEC_GRID}"

#-- input dir
INPUTDIR1=/data/lihl/code/3d-grid-generate/elliptic/project/output
#INPUTDIR1=/data/lihl/code/3d-grid-generate/parabolic/project/output
STRETCH_FILE1=`pwd`/arc_len_file1.txt
#STRETCH_FILE2=`pwd`/arc_len_file2.txt
#-- output and conf
PROJDIR=`pwd`/../project1
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
     "number_of_grid_points" : [541,541,301],
     "number_of_mpiprocs_in" : [2,2,2],
     "flag_stretch" : 0,
     "#stretch_file" : "${STRETCH_FILE1}"
   }
  ],
    
  "stretch_direction" : "z",
  "merge_direction" : "z",

  "number_of_mpiprocs_out" : [2,2,1],

  "check_orth" : 1,
  "check_jac" : 0,
  "check_step_xi" : 0,
  "check_step_et" : 0,
  "check_step_zt" : 0,
  "check_smooth_xi" : 1,
  "check_smooth_et" : 1,
  "check_smooth_zt" : 1,

  "flag_sample" : 0,
  "sample_factor_xi" : 2,
  "sample_factor_et" : 2,
  "sample_factor_zt" : 2,

  "grid_export_dir" : "${OUTPUT_DIR}",

  "flag_pml" : 1,
  "pml_layers" : {
         "number_of_pml_x1" : 20,
         "number_of_pml_x2" : 20,
         "number_of_pml_y1" : 20,
         "number_of_pml_y2" : 20,
         "number_of_pml_z1" : 20,
         "number_of_pml_z2" : 0
  }

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
