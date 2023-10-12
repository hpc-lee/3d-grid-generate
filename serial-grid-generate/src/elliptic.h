#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "gd_t.h"
#include "par_t.h"

int
diri_gene(gd_t *gdcurv, par_t *par);

int
set_src_diri(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, 
             float *P, float *Q, float *R, float *p_x, float *p_y, 
             float *p_z, float *g11_x, float *g22_y, float *g33_z, 
             float coef, int first_dire_itype, int second_dire_itype);

int
ghost_cal(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, 
          float *p_x, float *p_y, float *p_z,
          float *g11_x, float *g22_y, float *g33_z);

int
higen_gene(gd_t *gdcurv, par_t *par);

int
set_src_higen(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz,
              float *P, float *Q, float *R, float dx1, float dx2, 
              float dy1, float dy2, float dz1, float dz2,
              float coef, int first_dire_itype, int second_dire_itype);

int
interp_inner_source(float *P, float *Q, float *R,
                    int nx, int ny, int nz, float coef,
                    int first_dire_itype,int second_dire_itype);

#endif
