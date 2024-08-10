#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "constants.h"

#define PAR_MAX_STRLEN 1000
#define X_DIRE 1
#define Y_DIRE 2
#define Z_DIRE 3



typedef struct{

  int *num_of_points;
  int *num_of_procs_in;
  int *num_of_procs_out;
  int num_of_grid;
  char **import_dir;

  int *flag_stretch;
  char **stretch_file;

  char stretch_dire[PAR_MAX_STRLEN];
  int stretch_idire;

  char merge_dire[PAR_MAX_STRLEN];
  int merge_idire;

  int flag_pml;

  int  num_of_pml_x1;
  int  num_of_pml_x2;
  int  num_of_pml_y1;
  int  num_of_pml_y2;
  int  num_of_pml_z1;
  int  num_of_pml_z2;

  int grid_check;
  int check_orth;
  int check_jac;
  int check_step_xi;
  int check_step_et;
  int check_step_zt;
  int check_smooth_xi;
  int check_smooth_et;
  int check_smooth_zt;


  int flag_sample;
  int sample_factor_xi;
  int sample_factor_et;
  int sample_factor_zt;

  char export_dir[PAR_MAX_STRLEN];
  
} par_t;

int
par_read_from_file(char *par_fname, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

