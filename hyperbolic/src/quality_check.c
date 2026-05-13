#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "quality_check.h"
#include "constants.h"

int 
grid_quality_check(io_quality_t *io_quality, gd_t *gdcurv, par_t *par)
{
  if(par->check_orth == 1)
  {
    char quality_name1[100] = "orth_xiet";
    cal_xiet(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name1);
    char quality_name2[100] = "orth_xizt";
    cal_xizt(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name2);
    char quality_name3[100] = "orth_etzt";
    cal_etzt(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name3);
  }
  if(par->check_jac == 1)
  {
    char quality_name[100] = "jacobi";
    cal_jacobi(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }
  if(par->check_step_xi == 1)
  {
    char quality_name[100] = "step_xi";
    cal_step_xi(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }
  if(par->check_step_et == 1)
  {
    char quality_name[100] = "step_et";
    cal_step_et(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }
  if(par->check_step_zt == 1)
  {
    char quality_name[100] = "step_zt";
    cal_step_zt(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }
  if(par->check_smooth_xi == 1)
  {
    char quality_name[100] = "smooth_xi";
    cal_smooth_xi(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }
  if(par->check_smooth_et == 1)
  {
    char quality_name[100] = "smooth_et";
    cal_smooth_et(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }
  if(par->check_smooth_zt == 1)
  {
    char quality_name[100] = "smooth_zt";
    cal_smooth_zt(io_quality, gdcurv);
    quality_export(io_quality,par,quality_name);
  }

  return 0;
}

int 
cal_xiet(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;

  float trans = 180/PI;   // arc to angle
  float dot, len_xi, len_et, cos_angle;
  float x_xi, y_xi, z_xi, x_et, y_et, z_et;
  float var_min = 90;

  for(int k=0; k<nz-1; k++) {
    for(int j=0; j<ny-1; j++) {
      for(int i=0; i<nx-1; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + j*siz_iy + (i+1);
        // r_xi
        x_xi = x3d[iptr1] - x3d[iptr]; 
        y_xi = y3d[iptr1] - y3d[iptr]; 
        z_xi = z3d[iptr1] - z3d[iptr]; 

        iptr2 = k*siz_iz + (j+1)*siz_iy + i;
        // r_et
        x_et = x3d[iptr2] - x3d[iptr];
        y_et = y3d[iptr2] - y3d[iptr];
        z_et = z3d[iptr2] - z3d[iptr];
        
        // cos(theta) = a.b/(|a|*|b|)
        dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        cos_angle = dot/(len_xi*len_et); 
        // offset relative 90 degree
        var[iptr] = 90 - fabs((acos(cos_angle) * trans - 90));
        var_min = var_min < var[iptr] ? var_min : var[iptr];
      }
    }
  }
  fprintf(stdout,"xiet angle min is %f\n", var_min);
  fflush(stdout);
  
  extend_var(var, nx, ny, nz, siz_iy, siz_iz);

  return 0;
}

int 
cal_xizt(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;

  float trans = 180/PI;   // arc to angle
  float dot, len_xi, len_zt, cos_angle;
  float x_xi, y_xi, z_xi, x_zt, y_zt, z_zt;
  float var_min = 90;

  for(int k=0; k<nz-1; k++) {
    for(int j=0; j<ny-1; j++) {
      for(int i=0; i<nx-1; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + j*siz_iy + (i+1);
        // r_xi
        x_xi = x3d[iptr1] - x3d[iptr]; 
        y_xi = y3d[iptr1] - y3d[iptr]; 
        z_xi = z3d[iptr1] - z3d[iptr]; 

        iptr2 = (k+1)*siz_iz + j*siz_iy + i;
        // r_zt
        x_zt = x3d[iptr2] - x3d[iptr];
        y_zt = y3d[iptr2] - y3d[iptr];
        z_zt = z3d[iptr2] - z3d[iptr];
        
        // cos(theta) = a.b/(|a|*|b|)
        dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        cos_angle = dot/(len_xi*len_zt); 
        // offset relative 90 degree
        var[iptr] = 90 - fabs((acos(cos_angle) * trans - 90));
        var_min = var_min < var[iptr] ? var_min : var[iptr];
      }
    }
  }
  fprintf(stdout,"xizt angle min is %f\n", var_min);
  fflush(stdout);
  
  extend_var(var, nx, ny, nz, siz_iy, siz_iz);

  return 0;
}

int 
cal_etzt(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;

  float trans = 180/PI;   // arc to angle
  float dot, len_et, len_zt, cos_angle;
  float x_et, y_et, z_et, x_zt, y_zt, z_zt;
  float var_min = 90;

  for(int k=0; k<nz-1; k++) {
    for(int j=0; j<ny-1; j++) {
      for(int i=0; i<nx-1; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + (j+1)*siz_iy + i;
        // r_et
        x_et = x3d[iptr1] - x3d[iptr]; 
        y_et = y3d[iptr1] - y3d[iptr]; 
        z_et = z3d[iptr1] - z3d[iptr]; 

        iptr2 = (k+1)*siz_iz + j*siz_iy + i;
        // r_zt
        x_zt = x3d[iptr2] - x3d[iptr];
        y_zt = y3d[iptr2] - y3d[iptr];
        z_zt = z3d[iptr2] - z3d[iptr];
        
        // cos(theta) = a.b/(|a|*|b|)
        dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        cos_angle = dot/(len_et*len_zt); 
        // offset relative 90 degree
        var[iptr] = 90 - fabs((acos(cos_angle) * trans - 90));
        var_min = var_min < var[iptr] ? var_min : var[iptr];
      }
    }
  }
  fprintf(stdout,"etzt angle min is %f\n", var_min);
  fflush(stdout);
  
  
  extend_var(var, nx, ny, nz, siz_iy, siz_iz);

  return 0;
}

int 
cal_jacobi(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2,iptr3;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_xi, y_xi, z_xi;
  float x_et, y_et, z_et;
  float x_zt, y_zt, z_zt;
  float cros_x, cros_y, cros_z;
  // jacobi(volume) = (A X B).C 
  for(int k=0; k<nz-1; k++) {
    for(int j=0; j<ny-1; j++) {
      for(int i=0; i<nx-1; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + j*siz_iy + (i+1);
        // r_xi
        x_xi = x3d[iptr1] - x3d[iptr]; 
        y_xi = y3d[iptr1] - y3d[iptr]; 
        z_xi = z3d[iptr1] - z3d[iptr]; 

        iptr2 = k*siz_iz + (j+1)*siz_iy  + i;
        // r_et
        x_et = x3d[iptr2] - x3d[iptr];
        y_et = y3d[iptr2] - y3d[iptr];
        z_et = z3d[iptr2] - z3d[iptr];

        iptr3 = (k+1)*siz_iz + j*siz_iy  + i;
        // r_zt
        x_zt = x3d[iptr3] - x3d[iptr];
        y_zt = y3d[iptr3] - y3d[iptr];
        z_zt = z3d[iptr3] - z3d[iptr];

        // r_xi X r_et
        cros_x = y_xi*z_et - z_xi*y_et;
        cros_y = z_xi*x_et - x_xi*z_et;
        cros_z = x_xi*y_et - y_xi*x_et;

        var[iptr] = cros_x*x_zt + cros_y*y_zt + cros_z*z_zt;
      }
    }
  }

  extend_var(var, nx, ny, nz, siz_iy, siz_iz);

  return 0;
}

int 
cal_step_xi(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_xi, y_xi, z_xi;
  float len_xi;

  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx-1; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + j*siz_iy + (i+1);
        // r_xi
        x_xi = x3d[iptr1] - x3d[iptr]; 
        y_xi = y3d[iptr1] - y3d[iptr]; 
        z_xi = z3d[iptr1] - z3d[iptr]; 

        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        
        var[iptr] = len_xi;
      }
    }
  }

  // i = nx-1
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++)
    {
      int i = nx-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + j*siz_iy + (i-1);
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}

int 
cal_step_et(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_et, y_et, z_et;
  float len_et;

  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny-1; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + (j+1)*siz_iy + i;
        // r_et
        x_et = x3d[iptr1] - x3d[iptr]; 
        y_et = y3d[iptr1] - y3d[iptr]; 
        z_et = z3d[iptr1] - z3d[iptr]; 

        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        
        var[iptr] = len_et;
      }
    }
  }

  // j = ny-1
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      int j = ny-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + (j-1)*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}

int 
cal_step_zt(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_zt, y_zt, z_zt;
  float len_zt;

  for(int k=0; k<nz-1; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = (k+1)*siz_iz + j*siz_iy + i;
        // r_zt
        x_zt = x3d[iptr1] - x3d[iptr]; 
        y_zt = y3d[iptr1] - y3d[iptr]; 
        z_zt = z3d[iptr1] - z3d[iptr]; 

        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        
        var[iptr] = len_zt;
      }
    }
  }

  // k = nz-1
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++)
    {
      int k = nz-1;
      iptr  =  k*siz_iz + j*siz_iy + i;
      iptr1 = (k-1)*siz_iz + j*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}

int 
cal_smooth_xi(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_xi1, y_xi1, z_xi1;
  float x_xi2, y_xi2, z_xi2;
  float r1, r2, len_xi1, len_xi2;

  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + j*siz_iy + (i+1);
        iptr2 = k*siz_iz + j*siz_iy + (i-1);
        // r_xi
        x_xi1 = x3d[iptr1] - x3d[iptr]; 
        y_xi1 = y3d[iptr1] - y3d[iptr]; 
        z_xi1 = z3d[iptr1] - z3d[iptr]; 

        x_xi2 = x3d[iptr] - x3d[iptr2]; 
        y_xi2 = y3d[iptr] - y3d[iptr2]; 
        z_xi2 = z3d[iptr] - z3d[iptr2]; 

        len_xi1 = sqrt(pow(x_xi1,2) + pow(y_xi1,2) + pow(z_xi1,2));
        len_xi2 = sqrt(pow(x_xi2,2) + pow(y_xi2,2) + pow(z_xi2,2));

        r1 = len_xi1/len_xi2; 
        r2 = len_xi2/len_xi1; 
        
        var[iptr] = fmax(r1,r2);
      }
    }
  }

  // i = nx-1
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++)
    {
      int i = nx-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + j*siz_iy + (i-1);
      var[iptr] = var[iptr1];
    }
  }
  // i = 0
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++)
    {
      int i = 0;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + j*siz_iy + i+1;
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}

int 
cal_smooth_et(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_et1, y_et1, z_et1;
  float x_et2, y_et2, z_et2;
  float r1, r2, len_et1, len_et2;

  for(int k=0; k<nz; k++) {
    for(int j=1; j<ny-1; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = k*siz_iz + (j+1)*siz_iy + i;
        iptr2 = k*siz_iz + (j-1)*siz_iy + i;
        // r_et
        x_et1 = x3d[iptr1] - x3d[iptr]; 
        y_et1 = y3d[iptr1] - y3d[iptr]; 
        z_et1 = z3d[iptr1] - z3d[iptr]; 

        x_et2 = x3d[iptr] - x3d[iptr2]; 
        y_et2 = y3d[iptr] - y3d[iptr2]; 
        z_et2 = z3d[iptr] - z3d[iptr2]; 

        len_et1 = sqrt(pow(x_et1,2) + pow(y_et1,2) + pow(z_et1,2));
        len_et2 = sqrt(pow(x_et2,2) + pow(y_et2,2) + pow(z_et2,2));

        r1 = len_et1/len_et2; 
        r2 = len_et2/len_et1; 
        
        var[iptr] = fmax(r1,r2);
      }
    }
  }

  // j = ny-1
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      int j = ny-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + (j-1)*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }
  // j = 0
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      int j = 0;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + (j+1)*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}

int 
cal_smooth_zt(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *var = io_quality->var;
  float x_zt1, y_zt1, z_zt1;
  float x_zt2, y_zt2, z_zt2;
  float r1, r2, len_zt1, len_zt2;

  for(int k=1; k<nz-1; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = (k+1)*siz_iz + j*siz_iy + i;
        iptr2 = (k-1)*siz_iz + j*siz_iy + i;
        // r_zt
        x_zt1 = x3d[iptr1] - x3d[iptr]; 
        y_zt1 = y3d[iptr1] - y3d[iptr]; 
        z_zt1 = z3d[iptr1] - z3d[iptr]; 

        x_zt2 = x3d[iptr] - x3d[iptr2]; 
        y_zt2 = y3d[iptr] - y3d[iptr2]; 
        z_zt2 = z3d[iptr] - z3d[iptr2]; 

        len_zt1 = sqrt(pow(x_zt1,2) + pow(y_zt1,2) + pow(z_zt1,2));
        len_zt2 = sqrt(pow(x_zt2,2) + pow(y_zt2,2) + pow(z_zt2,2));

        r1 = len_zt1/len_zt2; 
        r2 = len_zt2/len_zt1; 
        
        var[iptr] = fmax(r1,r2);
      }
    }
  }

  // k = nz-1
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++)
    {
      int k = nz-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = (k-1)*siz_iz + j*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }
  // k = 0
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++)
    {
      int k = 0;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = (k+1)*siz_iz + j*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}

int extend_var(float *var, int nx, int ny, int nz,
               size_t siz_iy, size_t siz_iz)
{
  size_t iptr, iptr1;
  // i = nx-1
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++)
    {
      int i = nx-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + j*siz_iy + (i-1);
      var[iptr] = var[iptr1];
    }
  }

  // j = ny-1
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      int j = ny-1;
      iptr  = k*siz_iz + j*siz_iy + i;
      iptr1 = k*siz_iz + (j-1)*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }

  // k = nz-1
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++)
    {
      int k = nz-1;
      iptr  =  k*siz_iz + j*siz_iy + i;
      iptr1 = (k-1)*siz_iz + j*siz_iy + i;
      var[iptr] = var[iptr1];
    }
  }

  return 0;
}
