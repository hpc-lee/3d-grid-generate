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
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, int coef_x, int coef_y, int coef_z);

int
gd_info_set(gd_t *gdcurv, par_t *par, int iprocx, int iprocy, int iprocz,
           int *global_index, int *count);

int
cal_min_dist(gd_t *gdcurv, int *indx_i, int *indx_j, int *indx_k, float *dL_min);

#endif
