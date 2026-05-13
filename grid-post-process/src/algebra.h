#ifndef ALGEBRA_H
#define ALGEBRA_H


#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/
int 
xi_arc_stretch(gd_t *gdcurv, float *arc_len);

int 
et_arc_stretch(gd_t *gdcurv, float *arc_len);

int 
zt_arc_stretch(gd_t *gdcurv, float *arc_len);

int 
sample_interp(gd_t *gdcurv_new, gd_t *gdcurv);

float 
single_exp(float coef, float zeta);

#endif
