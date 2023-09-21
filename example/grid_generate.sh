#!/bin/bash

#set -x
set -e

date

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

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > ${PAR_FILE}
{

  "check_orth" : 1,
  "check_jac" : 1,
  "check_step_x" : 1,
  "check_step_y" : 1,
  "check_step_z" : 1,
  "check_smooth_x" : 1,
  "check_smooth_y" : 1,
  "check_smooth_z" : 1,

  "flag_strech_x" : 0,
  "strech_x_coef" : 0.0001,
  "flag_strech_y" : 0,
  "strech_y_coef" : 0.0001,
  "flag_strech_z" : 0,
  "strech_z_coef" : 0.0001,

  "flag_sample_x" : 0,
  "sample_factor_x" : 2.0,
  "flag_sample_y" : 0,
  "sample_factor_y" : 2.0,
  "flag_sample_z" : 0,
  "sample_factor_z" : 1.0,

  "geometry_input_file" : "${INPUTDIR}/data_file_3d.txt",
  "grid_export_dir" : "${OUTPUT_DIR}",

  "grid_method" : {
      "#linear_TFI" : "",
      "#elli_diri" : {
          "coef" : -20,
          "iter_err" : 1E-2,
          "max_iter" : 5E3,
          "first_dire"  : "z",
          "second_dire" : "y"
      },
      "#elli_higen" : {
          "coef" : -20,
          "iter_err" : 1E-2,
          "max_iter" : 5E3,
          "distance" : [100,100,100,100,100,100],
          "first_dire"  : "z",
          "second_dire" : "y"
      },
      "#parabolic" : {
          "coef" : -50,
          "direction" : "z",
          "o2i" : 1
      },
      "hyperbolic" : {
          "coef" : 20,
          "bdry_x_type" : 2,
          "epsilon_x" : 0,
          "bdry_y_type" : 2,
          "epsilon_y" : 0,
          "direction" : "z",
          "o2i" : 1,
          "step_input_file" : "${INPUTDIR}/step_file_3d.txt"
      }
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
