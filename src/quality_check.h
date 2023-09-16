#ifndef QUALITY_CHECK_H
#define QUALITY_CHECK_H

#include "gd_t.h"
#include "par_t.h"
#include "io_funcs.h"

int 
grid_quality_check(io_quality_t *io_quality, gd_t *gdcurv, par_t *par);

int 
cal_xiet(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_xizt(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_etzt(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_jacobi(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_x(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_y(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_z(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_x(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_y(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_z(io_quality_t *io_quality, gd_t *gdcurv);

int extend_var(float *var, int nx, int ny, int nz,
               size_t siz_iy, size_t siz_iz);

#endif
