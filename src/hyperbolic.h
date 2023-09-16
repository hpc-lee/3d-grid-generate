#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H

#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
hyper_gene(gd_t *gdcurv, float coef, int o2i, int bdry_itype, float epsilon);

int
cal_smooth_coef(float coef, float *x3d, float *z3d, int nx, int nz, int k, float *coef_e);

int 
cal_matrix(float *x3d, float *z3d, int nx, int k, float *step,
           double *a, double *b, double *c, double *d, float *area);

int
modify_smooth(float *x3d, float *z3d, int nx, int k, double *a,
              double *b, double *c, double *d, float *coef_e);

int
modify_bdry(int n, double *a, double *b, double *c, double *d,
            float epsilon, int bdry_itype);

int
assign_coords(double *xz, float *x3d, float *z3d, int nx, int k,
              float epsilon, int bdry_itype);

 #endif
