#ifndef GD_CURV_H
#define GD_CURV_H

#include "par_t.h"
#include "mympi_t.h"
/*************************************************
 * structure
 *************************************************/

typedef struct {

  int ni, nj, nk;
  int nx, ny, nz;
  int ni1, ni2;
  int nj1, nj2;
  int nk1, nk2;
  // global index
  int gni1, gnj1, gnk1; // global index, do not accout ghost point
  int gni2, gnj2, gnk2; // global index

  int total_nx;
  int total_ny;
  int total_nz;

  int ncmp;
  
  float *v4d; // pointer to var
  float *x3d; 
  float *y3d; 
  float *z3d;

  float *v4d_tmp; // temp space pointer to var
  float *x3d_tmp; 
  float *y3d_tmp; 
  float *z3d_tmp;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  char fname_part[CONST_MAX_STRLEN];
  char output_dir[CONST_MAX_STRLEN];

} gd_t;

typedef struct {

  float *var;
  float *x1;
  float *x2;
  float *y1;
  float *y2;
  float *z1;
  float *z2;
  int total_nx;
  int total_ny;
  int total_nz;

} bdry_t;
/*************************************************
 * function prototype
 *************************************************/
int
gd_info_set(gd_t *gdcurv, mympi_t *mympi,
            par_t *par);

int
set_output_dir(gd_t *gdcurv,
               mympi_t *mympi,
               char *output_dir);

int 
init_gdcurv(gd_t *gdcurv);

int
init_bdry(bdry_t *bdry, par_t *par);

int
read_bdry(int myid, bdry_t *bdry, char *geometry_file);

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z);

int 
check_bdry(float *x1, float *x2, float *y1, float *y2, 
           float *z1, float *z2, int nx, int ny, int nz);

int
gd_info_print(gd_t *gdcurv, mympi_t *mympi);

int
gd_curv_coord_exchange(gd_t *gdcurv, int *neighid, MPI_Comm topocomm);

int
grid_mesg_init(mympi_t *mympi, gd_t *gdcurv);

int
grid_pack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x3d, float *y3d, float *z3d);

int
grid_unpack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x3d, float *y3d, float *z3d);

int
cal_min_dist(gd_t *gdcurv, int *indx_i, int *indx_j, int *indx_k, float *dL_min);

#endif
