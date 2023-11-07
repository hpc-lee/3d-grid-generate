#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "constants.h"

#define PAR_MAX_STRLEN 1000


typedef struct{

  int number_of_grid_points_x;
  int number_of_grid_points_y;
  int number_of_grid_points_z;

  int number_of_mpiprocs_x_in;
  int number_of_mpiprocs_y_in;
  int number_of_mpiprocs_z_in;

  int number_of_mpiprocs_x_out;
  int number_of_mpiprocs_y_out;
  int number_of_mpiprocs_z_out;

  int  number_of_pml_x1;
  int  number_of_pml_x2;
  int  number_of_pml_y1;
  int  number_of_pml_y2;
  int  number_of_pml_z1;
  int  number_of_pml_z2;

  int grid_check;
  int check_orth;
  int check_jac;
  int check_step_xi;
  int check_step_et;
  int check_step_zt;
  int check_smooth_xi;
  int check_smooth_et;
  int check_smooth_zt;

  int flag_strech_xi;
  int flag_strech_et;
  int flag_strech_zt;
  float strech_xi_coef;
  float strech_et_coef;
  float strech_zt_coef;

  int flag_sample;
  int sample_factor_xi;
  int sample_factor_et;
  int sample_factor_zt;

  char import_dir[PAR_MAX_STRLEN];
  char export_dir[PAR_MAX_STRLEN];
  
} par_t;

int
par_read_from_file(char *par_fname, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

