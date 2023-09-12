#ifndef GD_CURV_H
#define GD_CURV_H

/*************************************************
 * structure
 *************************************************/

typedef struct {

  int nx;
  int nz;
  int ncmp;
  
  float *v3d; // pointer to var
  float *x2d; 
  float *z2d;
  
  float *step; // for hyperbolic

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

} gd_t;

/*************************************************
 * function prototype
 *************************************************/

int 
init_gdcurv(gd_t *gdcurv, int nx, int nz);

int 
grid_init_set(gd_t *gdcurv, char *input_file);

int
grid_init_set_hyper(gd_t *gdcurv, char *geometry_file, char *step_file);

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z);

int 
check_bdry(float *x1, float *x2, float *z1, float *z2, int nx, int nz);

int
flip_coord(float *coord, int nx, int nz);

int
permute_coord(gd_t *gdcurv);

#endif
