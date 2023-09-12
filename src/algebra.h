#ifndef ALGEBRA_H
#define ALGEBRA_H


#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
linear_tfi(gd_t *gdcurv);

int
one_hermite(gd_t *gdcurv, float coef);

int 
zt_arc_strech(gd_t *gdcurv, float coef);

int 
xi_arc_strech(gd_t *gdcurv, float coef);

int 
sample_interp(gd_t *gdcurv_new, gd_t *gdcurv);

float 
single_exp(float coef, float zeta);

#endif
