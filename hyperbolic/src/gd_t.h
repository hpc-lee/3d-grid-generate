#ifndef GD_CURV_H
#define GD_CURV_H

#include "par_t.h"
/*************************************************
 * structure
 *************************************************/

typedef struct {

  int nx;
  int ny;
  int nz;
  int ncmp;
  
  float *v4d; // pointer to var
  float *x3d; 
  float *y3d; 
  float *z3d;
  
  float *step; // for hyperbolic

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;
} gd_t;

/*************************************************
 * function prototype
 *************************************************/

int 
init_gdcurv(gd_t *gdcurv, int nx, int ny, int nz);

int
grid_init_set_hyper(gd_t *gdcurv, par_t *par);

int
flip_coord_z(gd_t *gdcurv);

int
permute_coord_x(gd_t *gdcurv);

int
permute_coord_y(gd_t *gdcurv);

int
cal_min_dist(gd_t *gdcurv, int *indx_i, int *indx_j, int *indx_k, float *dL_min);

#endif
