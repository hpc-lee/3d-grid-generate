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
    char quality_name[100] = "orth";
    cal_orth(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }
  if(par->check_jac == 1)
  {
    char quality_name[100] = "jacobi";
    cal_jacobi(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }
  if(par->check_ratio == 1)
  {
    char quality_name[100] = "ratio";
    cal_ratio(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }
  if(par->check_step_x == 1)
  {
    char quality_name[100] = "step_x";
    cal_step_x(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }
  if(par->check_step_z == 1)
  {
    char quality_name[100] = "step_z";
    cal_step_z(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }
  if(par->check_smooth_x == 1)
  {
    char quality_name[100] = "smooth_x";
    cal_smooth_x(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }
  if(par->check_smooth_z == 1)
  {
    char quality_name[100] = "smooth_z";
    cal_smooth_z(io_quality, gdcurv);
    quality_export(io_quality,par->grid_export_dir,quality_name);
  }

  return 0;
}

int 
cal_orth(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;

  float trans = 180/PI;   // arc to angle
  float dot, len_xi, len_zt, cos_angle;
  float x_xi, z_xi, x_zt, z_zt;

  for(int k=0; k<nz-1; k++) {
    for(int i=0; i<nx-1; i++)
    {
      iptr  = k*nx + i;
      iptr1 = k*nx + (i+1);
      // r_xi
      x_xi = x2d[iptr1] - x2d[iptr]; 
      z_xi = z2d[iptr1] - z2d[iptr]; 

      iptr2 = (k+1)*nx + i;
      // r_zt
      x_zt = x2d[iptr2] - x2d[iptr];
      z_zt = z2d[iptr2] - z2d[iptr];
      
      // cos(theta) = a.b/(|a|*|b|)
      dot = x_xi*x_zt + z_xi*z_zt;
      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
      cos_angle = dot/(len_xi*len_zt); 
      // offset relative 90 degree
      var[iptr] = fabs((acos(cos_angle) * trans - 90));
    }
  }
  
  // i = nx-1
  for(int k=0; k<nz; k++)
  {
    int i = nx-1;
    iptr  = k*nx + i;
    iptr1 = k*nx + (i-1);
    var[iptr] = var[iptr1];
  }

  // k = nz-1
  for(int i=0; i<nx; i++)
  {
    int k = nz-1;
    iptr  =  k*nx + i;
    iptr1 = (k-1)*nx + i;
    var[iptr] = var[iptr1];
  }

  return 0;
}

int 
cal_jacobi(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;
  float x_xi, z_xi, x_zt, z_zt;
  
  for(int k=0; k<nz-1; k++) {
    for(int i=0; i<nx-1; i++)
    {
      iptr  = k*nx + i;
      iptr1 = k*nx + (i+1);
      // r_xi
      x_xi = x2d[iptr1] - x2d[iptr]; 
      z_xi = z2d[iptr1] - z2d[iptr]; 

      iptr2 = (k+1)*nx + i;
      // r_zt
      x_zt = x2d[iptr2] - x2d[iptr];
      z_zt = z2d[iptr2] - z2d[iptr];

      var[iptr] = x_xi*z_zt - z_xi*x_zt;
    }
  }

  // i = nx-1
  for(int k=0; k<nz; k++)
  {
    int i = nx-1;
    iptr  = k*nx + i;
    iptr1 = k*nx + (i-1);
    var[iptr] = var[iptr1];
  }

  // k = nz-1
  for(int i=0; i<nx; i++)
  {
    int k = nz-1;
    iptr  =  k*nx + i;
    iptr1 = (k-1)*nx + i;
    var[iptr] = var[iptr1];
  }

  return 0;
}

int 
cal_ratio(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;
  float x_xi, z_xi, x_zt, z_zt;
  float r1, r2, len_xi, len_zt;

  for(int k=0; k<nz-1; k++) {
    for(int i=0; i<nx-1; i++)
    {
      iptr  = k*nx + i;
      iptr1 = k*nx + (i+1);
      // r_xi
      x_xi = x2d[iptr1] - x2d[iptr]; 
      z_xi = z2d[iptr1] - z2d[iptr]; 

      iptr2 = (k+1)*nx + i;
      // r_zt
      x_zt = x2d[iptr2] - x2d[iptr];
      z_zt = z2d[iptr2] - z2d[iptr];

      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));

      r1 = len_xi/len_zt; 
      r2 = len_zt/len_xi; 
      
      var[iptr] = fmax(r1,r2);
    }
  }

  // i = nx-1
  for(int k=0; k<nz; k++)
  {
    int i = nx-1;
    iptr  = k*nx + i;
    iptr1 = k*nx + (i-1);
    var[iptr] = var[iptr1];
  }

  // k = nz-1
  for(int i=0; i<nx; i++)
  {
    int k = nz-1;
    iptr  =  k*nx + i;
    iptr1 = (k-1)*nx + i;
    var[iptr] = var[iptr1];
  }

  return 0;
}

int 
cal_step_x(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;
  float x_xi, z_xi;
  float len_xi;

  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx-1; i++)
    {
      iptr  = k*nx + i;
      iptr1 = k*nx + (i+1);
      // r_xi
      x_xi = x2d[iptr1] - x2d[iptr]; 
      z_xi = z2d[iptr1] - z2d[iptr]; 

      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      
      var[iptr] = len_xi;
    }
  }

  // i = nx-1
  for(int k=0; k<nz; k++)
  {
    int i = nx-1;
    iptr  = k*nx + i;
    iptr1 = k*nx + (i-1);
    var[iptr] = var[iptr1];
  }

  return 0;
}

int 
cal_step_z(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;
  float x_zt, z_zt;
  float len_zt;

  for(int k=0; k<nz-1; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr  = k*nx + i;
      iptr1 = (k+1)*nx + i;
      // r_xi
      x_zt = x2d[iptr1] - x2d[iptr]; 
      z_zt = z2d[iptr1] - z2d[iptr]; 

      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
      
      var[iptr] = len_zt;
    }
  }

  // k = nz-1
  for(int i=0; i<nx; i++)
  {
    int k = nz-1;
    iptr  = k*nx + i;
    iptr1 = (k-1)*nx + i;
    var[iptr] = var[iptr1];
  }

  return 0;
}

int 
cal_smooth_x(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;
  float x_xi1, z_xi1, x_xi2, z_xi2;
  float r1, r2, len_xi1, len_xi2;

  for(int k=0; k<nz; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = k*nx + i;
      iptr1 = k*nx + (i-1);
      iptr2 = k*nx + (i+1);
      // r_xi
      x_xi1 = x2d[iptr] - x2d[iptr1]; 
      z_xi1 = z2d[iptr] - z2d[iptr1]; 

      x_xi2 = x2d[iptr2] - x2d[iptr]; 
      z_xi2 = z2d[iptr2] - z2d[iptr]; 

      len_xi1 = sqrt(pow(x_xi1,2) + pow(z_xi1,2));
      len_xi2 = sqrt(pow(x_xi2,2) + pow(z_xi2,2));

      r1 = len_xi1/len_xi2; 
      r2 = len_xi2/len_xi1; 
      
      var[iptr] = fmax(r1,r2);
    }
  }

  // i = 0
  for(int k=0; k<nz; k++)
  {
    int i = 0;
    iptr  = k*nx + i;
    iptr1 = k*nx + (i+1);
    var[iptr] = var[iptr1];
  }
  // i = nx-1
  for(int k=0; k<nz; k++)
  {
    int i = nx-1;
    iptr  = k*nx + i;
    iptr1 = k*nx + (i-1);
    var[iptr] = var[iptr1];
  }

  return 0;
}

int 
cal_smooth_z(io_quality_t *io_quality, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr,iptr1,iptr2;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *var = io_quality->var;
  float x_zt1, z_zt1, x_zt2, z_zt2;
  float r1, r2, len_zt1, len_zt2;

  for(int k=1; k<nz-1; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr  = k*nx + i;
      iptr1 = (k-1)*nx + i;
      iptr2 = (k+1)*nx + i;
      // r_xi
      x_zt1 = x2d[iptr] - x2d[iptr1]; 
      z_zt1 = z2d[iptr] - z2d[iptr1]; 

      x_zt2 = x2d[iptr2] - x2d[iptr]; 
      z_zt2 = z2d[iptr2] - z2d[iptr]; 

      len_zt1 = sqrt(pow(x_zt1,2) + pow(z_zt1,2));
      len_zt2 = sqrt(pow(x_zt2,2) + pow(z_zt2,2));

      r1 = len_zt1/len_zt2; 
      r2 = len_zt2/len_zt1; 
      
      var[iptr] = fmax(r1,r2);
    }
  }

  // k = 0
  for(int i=0; i<nx; i++)
  {
    int k = 0;
    iptr  = k*nx + i;
    iptr1 = (k+1)*nx + i;
    var[iptr] = var[iptr1];
  }
  // k = nz-1
  for(int i=0; i<nx; i++)
  {
    int k = nz-1;
    iptr  = k*nx + i;
    iptr1 = (k-1)*nx + i;
    var[iptr] = var[iptr1];
  }

  return 0;
}
