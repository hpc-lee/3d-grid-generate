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
grid_init_set(gd_t *gdcurv, char *input_file);

int
grid_init_set_hyper(gd_t *gdcurv, par_t *par);

int 
check_bdry(float *x1, float *x2, float *y1, float *y2,float *z1, float *z2,
           int nx, int ny, int nz);

int
flip_coord_z(gd_t *gdcurv);

int
permute_coord_x(gd_t *gdcurv);

int
permute_coord_y(gd_t *gdcurv);

int
gd_info_set(par_t *par, int iprocx, int iprocy, int iprocz,
           int *global_index, int *count);

#endif
