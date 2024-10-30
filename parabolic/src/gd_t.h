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
  int gni1, gni2; // global index, do not accout ghost point
  int gnj1, gnj2; 
  int gnk1, gnk2; 

  int total_nx;
  int total_ny;
  int total_nz;

  int ncmp;
  
  float *v4d; // pointer to var
  float *x3d; 
  float *y3d; 
  float *z3d;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  float *step;

  char fname_part[CONST_MAX_STRLEN];
  char output_dir[CONST_MAX_STRLEN];

} gd_t;


/*************************************************
 * function prototype
 *************************************************/

int
init_gdcurv(gd_t *gdcurv);

int
gd_info_set(gd_t *gdcurv, mympi_t *mympi,
            par_t *par);

int
gd_info_print(gd_t *gdcurv, mympi_t *mympi);

int
set_output_dir(gd_t *gdcurv, mympi_t *mympi,
               par_t *par);

int
read_bdry_file(gd_t *gdcurv, par_t *par);

int
permute_coord_x(gd_t *gdcurv);

int
permute_coord_y(gd_t *gdcurv);

int
flip_coord_z(gd_t *gdcurv);

int
grid_mesg_init(mympi_t *mympi, gd_t *gdcurv);

int
grid_pack_mesg(mympi_t *mympi, gd_t *gdcurv, int k);

int
grid_unpack_mesg(mympi_t *mympi, gd_t *gdcurv, int k);

int
cal_min_dist(gd_t *gdcurv, int *indx_i, int *indx_j, int *indx_k, float *dL_min);

#endif
