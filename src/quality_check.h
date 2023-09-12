#ifndef QUALITY_CHECK_H
#define QUALITY_CHECK_H

#include "gd_t.h"
#include "par_t.h"
#include "io_funcs.h"

int 
grid_quality_check(io_quality_t *io_quality, gd_t *gdcurv, par_t *par);


int 
cal_orth(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_jacobi(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_ratio(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_x(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_z(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_x(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_z(io_quality_t *io_quality, gd_t *gdcurv);

#endif
