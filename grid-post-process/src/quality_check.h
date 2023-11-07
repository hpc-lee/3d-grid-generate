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
cal_step_xi(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_et(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_step_zt(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_xi(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_et(io_quality_t *io_quality, gd_t *gdcurv);

int 
cal_smooth_zt(io_quality_t *io_quality, gd_t *gdcurv);

int extend_var(float *var, int nx, int ny, int nz,
               size_t siz_iy, size_t siz_iz);

#endif
