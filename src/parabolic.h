#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
para_gene(gd_t *gdcurv,float coef, int o2i);

int 
predict_point(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, int k, int o2i, 
              float coef, float *x1_len, float *x2_len, float *y1_len, float *y2_len);

int
update_point(float *x3d, float *y3d, float *z3d, float *thomas, int nx, int ny, int nz, int k);

int
bdry_effct(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, int k);

int
cal_bdry_arc_length(float *x3d, float *y3d, float *z3d, int nx,
                    int ny, int nz, float *x1_len, float *x2_len,
                    float *y1_len, float *y2_len);

#endif
