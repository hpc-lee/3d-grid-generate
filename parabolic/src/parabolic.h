#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "gd_t.h"

/*************************************************
 * function prototype
 *************************************************/

int 
para_gene(gd_t *gdcurv, mympi_t *mympi, par_t *par);

int 
predict_point(gd_t *gdcurv, mympi_t *mympi,
              int k, int t2b, float coef, 
              float *step_len);

int
update_point(gd_t *gdcurv, float *var_th, int k, float *coord); 

int
exchange_coord(gd_t *gdcurv, mympi_t *mympi, int k, int num_of_s_reqs, int num_of_r_reqs);

int
flip_step_z(float *step, int nz);

#endif
