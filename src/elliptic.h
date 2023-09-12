#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "gd_t.h"
#include "par_t.h"

int
diri_gene(gd_t *gdcurv, par_t *par);

int
set_src_diri(float *x2d, float *z2d, int nx, int nz, 
             float *P, float *Q, float *p_x, float *p_z,
             float *g11_x, float *g22_z, float coef);

int
ghost_cal(float *x2d, float *z2d, int nx, int nz, float *p_x, float *p_z,
          float *g11_x, float *g22_z);

int
higen_gene(gd_t *gdcurv, par_t *par);

int
set_src_higen(float *x2d, float *z2d, int nx, int nz,
              float *P, float *Q, float coef,
              float dx1, float dx2, float dz1, float dz2);

#endif
