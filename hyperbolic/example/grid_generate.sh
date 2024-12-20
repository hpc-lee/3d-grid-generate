#!/bin/bash

#set -x
set -e

date

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

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{
  "number_of_grid_points_x" : 541,
  "number_of_grid_points_y" : 541,
  "number_of_grid_points_z" : 301,

  "check_orth" : 1,
  "check_jac" : 0,
  "check_step_xi" : 0,
  "check_step_et" : 0,
  "check_step_zt" : 0,
  "check_smooth_xi" : 1,
  "check_smooth_et" : 1,
  "check_smooth_zt" : 1,

  "geometry_input_file" : "${INPUTDIR}/data_file_3d.txt",
  "grid_export_dir" : "${OUTPUT_DIR}",

  "hyperbolic" : {
      "coef" : 40.0,
      "bdry_x_type" : 1,
      "epsilon_x" : 0.0,
      "bdry_y_type" : 1,
      "epsilon_y" : 0.0,
      "direction" : "z",
      "o2i" : 1,
      "step_input_file" : "${INPUTDIR}/step_file_3d.txt"
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

printf "\nStart grid generate ...\n";
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
