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

  size_t *cmp_pos;
  char  **cmp_name;

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
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_y, float coef_z);

int 
check_bdry(float *x1, float *x2, float *y1, float *y2,float *z1, float *z2,
           int nx, int ny, int nz);

int
flip_coord_z(gd_t *gdcurv);

int
permute_coord_x(gd_t *gdcurv);

int
permute_coord_y(gd_t *gdcurv);

#endif
