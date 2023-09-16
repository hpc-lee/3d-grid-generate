#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "constants.h"
#include "par_t.h"

#define PAR_MAX_STRLEN 1000
#define TFI 1
#define ELLI_DIRI 2
#define ELLI_HIGEN 3
#define PARABOLIC 4
#define HYPERBOLIC 5

#define X_DIRE 1
#define Y_DIRE 2
#define Z_DIRE 3

typedef struct{

  int grid_check;
  int check_orth;
  int check_jac;
  int check_step_x;
  int check_step_y;
  int check_step_z;
  int check_smooth_x;
  int check_smooth_y;
  int check_smooth_z;

  int flag_strech_x;
  int flag_strech_y;
  int flag_strech_z;
  float strech_x_coef;
  float strech_y_coef;
  float strech_z_coef;

  int flag_sample_x;
  int flag_sample_y;
  int flag_sample_z;
  float sample_factor_x;
  float sample_factor_y;
  float sample_factor_z;

  char geometry_input_file[PAR_MAX_STRLEN];
  char step_input_file[PAR_MAX_STRLEN];
  char grid_export_dir[PAR_MAX_STRLEN];
  
  // TFI elliptic-dirichlet
  // elliptic-hilgenstock
  // parabolic hyperbolic
  int method_itype;

  int first_dire_itype;   // first second for elliptic
  int second_dire_itype;
  char first_dire[PAR_MAX_STRLEN];
  char second_dire[PAR_MAX_STRLEN];

  int dire_itype;  // parabolic hyperbolic
  char direction[PAR_MAX_STRLEN];
  float coef;

  float distance[6];  // for higenstock dx1,dx2,dy1,dy2,dz1,dz2 
  float iter_err;   // iteration error
  int max_iter;  // max iterations
  int o2i;  // for parabolic hyperbolic. outer to inner
  int bdry_x_itype; // for hyperbolic
  int bdry_y_itype; // for hyperbolic
  float epsilon_x;  // for hyperbolic
  float epsilon_y;  // for hyperbolic
} par_t;

int
par_read_from_file(char *par_fname, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

