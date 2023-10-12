#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "cJSON.h"
#include "constants.h"

#define PAR_MAX_STRLEN 1000
#define ELLI_DIRI 1
#define ELLI_HIGEN 2

#define X_DIRE 1
#define Z_DIRE 2

typedef struct{

  //grid size
  int number_of_grid_points_x;
  int number_of_grid_points_y;
  int number_of_grid_points_z;

  // MPI
  int number_of_mpiprocs_x;
  int number_of_mpiprocs_y;
  int number_of_mpiprocs_z;

  // number of pml
  int number_of_pml_x1;
  int number_of_pml_x2;
  int number_of_pml_y1;
  int number_of_pml_y2;
  int number_of_pml_z1;
  int number_of_pml_z2;

  int grid_check;
  int check_orth;
  int check_jac;
  int check_step_xi;
  int check_step_et;
  int check_step_zt;
  int check_smooth_xi;
  int check_smooth_et;
  int check_smooth_zt;

  char geometry_input_file[PAR_MAX_STRLEN];
  char output_dir[PAR_MAX_STRLEN];
  
  // TFI hermite elliptic-dirichlet
  // elliptic-hilgenstock
  // parabolic hyperbolic
  int method_itype;

  float coef;

  float distance[6];  // for higenstock dx1,dx2,dy1,dy2,dz1,dz2 
  float iter_err;   // iteration error
  int max_iter;  // max iterations
} par_t;

int
par_mpi_get(char *par_fname, int myid, MPI_Comm comm, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

