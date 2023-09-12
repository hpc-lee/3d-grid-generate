#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "algebra.h"
#include "lib_mem.h"

// linear tfi interpolation
// U means x-direction
// W means z-direction
int linear_tfi(gd_t *gdcurv)
{
  int nx = gdcurv->nx; 
  int nz = gdcurv->nz; 
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float xi,zeta;
  float a0,a1,c0,c1;
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  size_t iptr5,iptr6,iptr7,iptr8;
  float U_x,W_x,UW_x;
  float U_z,W_z,UW_z;
  
  for (int k=1; k<nz-1; k++) {
    for (int i=1; i<nx-1; i++)
    {
      xi   = (1.0*i)/(nx-1); // 1.0* equal int to float
      zeta = (1.0*k)/(nz-1);
      a0 = 1-xi;
      a1 = xi;
      c0 = 1-zeta;
      c1 = zeta;

      iptr = k*nx + i;  //(i,k) 
      // 4 boudary points
      iptr1 = k*nx;     //(0,k)
      iptr2 = k*nx + nx-1; //(nx-1,k)
      iptr3 = i;         //(i,0)
      iptr4 = (nz-1)*nx + i; //(i,nz-1)
      // 4 corner points
      iptr5 = 0;   //(0,0)
      iptr6 = (nz-1)*nx;   //(0,nz-1)
      iptr7 = nx-1;   //(nx-1,0)
      iptr8 = (nz-1)*nx + nx-1;   //(nx-1,nz-1)

      U_x = a0*x2d[iptr1] + a1*x2d[iptr2];
      W_x = c0*x2d[iptr3] + c1*x2d[iptr4];
      UW_x = a0*c0*x2d[iptr5] + a0*c1*x2d[iptr6] \
           + a1*c0*x2d[iptr7] + a1*c1*x2d[iptr8];

      U_z = a0*z2d[iptr1] + a1*z2d[iptr2];
      W_z = c0*z2d[iptr3] + c1*z2d[iptr4];
      UW_z = a0*c0*z2d[iptr5] + a0*c1*z2d[iptr6] \
           + a1*c0*z2d[iptr7] + a1*c1*z2d[iptr8];

      x2d[iptr] = U_x + W_x - UW_x;
      z2d[iptr] = U_z + W_z - UW_z;
    }
  }
  
  return 0;
}

// hermite interpolation
// four boundary x1(left), x2(right), z1(bottom) ,z2(top)
int
one_hermite(gd_t *gdcurv, float coef)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float zeta,c0,c1,c10,c11;
  float tan_z1_x,tan_z1_z,tan_z2_x,tan_z2_z;
  float x_len,z_len;
  float abs_z1, abs_z2;
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float *dh_len;
  float *nor_z1_x;
  float *nor_z1_z;
  float *nor_z2_x;
  float *nor_z2_z;

  dh_len = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z1_x = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z1_z = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z2_x = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z2_z = (float *)mem_calloc_1d_float(nx, 0.0, "init");

  for (int i=0; i<nx; i++)
  {
    iptr1 = (nz-1)*nx + i; //(i,nz-1)
    iptr2 = i;             //(i,1)
    x_len = x2d[iptr1] - x2d[iptr2];
    z_len = z2d[iptr1] - z2d[iptr2];
    dh_len[i] = sqrt(pow(x_len,2) + pow(z_len,2));
  }

  for (int i=1; i<nx-1; i++) 
  {
    iptr1 = i+1;   //(i+1,1)
    iptr2 = i-1;   //(i-1,1)
    iptr3 = (nz-1)*nx + (i+1);   //(i+1,nz-1)
    iptr4 = (nz-1)*nx + (i-1);   //(i-1,nz-1)
    tan_z1_x = (x2d[iptr1] - x2d[iptr2])/2;
    tan_z1_z = (z2d[iptr1] - z2d[iptr2])/2;
    tan_z2_x = (x2d[iptr3] - x2d[iptr4])/2;
    tan_z2_z = (z2d[iptr3] - z2d[iptr4])/2;
    nor_z1_x[i] = -tan_z1_z;
    nor_z1_z[i] = tan_z1_x;
    nor_z2_x[i] = -tan_z2_z;
    nor_z2_z[i] = tan_z2_x;
  }
  
  // i=0 
  iptr1 = 1;   //(1,1)
  iptr2 = 0;   //(0,1)
  iptr3 = (nz-1)*nx + 1;   //(1,nz-1)
  iptr4 = (nz-1)*nx + 0;   //(0,nz-1)
  nor_z1_x[0] = -(z2d[iptr1]-z2d[iptr2]);
  nor_z1_z[0] = x2d[iptr1]-x2d[iptr2];
  nor_z2_x[0] = -(z2d[iptr3]-z2d[iptr4]);
  nor_z2_z[0] = x2d[iptr3]-x2d[iptr4];
  // i=nx-1
  iptr1 = nx-1;   //(nx-1,1)
  iptr2 = nx-2;   //(nx-2,1)
  iptr3 = (nz-1)*nx + nx-1;   //(nx-1,nz-1)
  iptr4 = (nz-1)*nx + nx-2;   //(nx-2,nz-1)
  nor_z1_x[nx-1] = -(z2d[iptr1]-z2d[iptr2]);
  nor_z1_z[nx-1] = x2d[iptr1]-x2d[iptr2];
  nor_z2_x[nx-1] = -(z2d[iptr3]-z2d[iptr4]);
  nor_z2_z[nx-1] = x2d[iptr3]-x2d[iptr4];

  // cal 1st order derivative coef, this effect 
  // boundary orth length
  for (int i=0; i<nx; i++) 
  {
    abs_z1 = sqrt(pow(nor_z1_x[i],2) + pow(nor_z1_z[i],2));
    abs_z2 = sqrt(pow(nor_z2_x[i],2) + pow(nor_z2_z[i],2));
    nor_z1_x[i] = coef*dh_len[i]*(nor_z1_x[i]/abs_z1);
    nor_z1_z[i] = coef*dh_len[i]*(nor_z1_z[i]/abs_z1);
    nor_z2_x[i] = coef*dh_len[i]*(nor_z2_x[i]/abs_z2);
    nor_z2_z[i] = coef*dh_len[i]*(nor_z2_z[i]/abs_z2);
  }

  for (int i=0; i<nx; i++) {
    for (int k=1; k<nz-1; k++)
    {
      zeta = (1.0*k)/(nz-1);
      c0 =  2*pow(zeta,3) - 3*pow(zeta,2) + 1;
      c1 = -2*pow(zeta,3) + 3*pow(zeta,2);
      c10 = pow(zeta,3)-2*pow(zeta,2)+zeta;
      c11 = pow(zeta,3)-pow(zeta,2);
      
      iptr = k*nx + i;
      iptr1 = i;
      iptr2 = (nz-1)*nx+i;
      x2d[iptr] = c0*x2d[iptr1] + c1*x2d[iptr2] + c10*nor_z1_x[i] + c11*nor_z2_x[i];
      z2d[iptr] = c0*z2d[iptr1] + c1*z2d[iptr2] + c10*nor_z1_z[i] + c11*nor_z2_z[i];
    }
  }

  free(dh_len);
  free(nor_z1_x); 
  free(nor_z1_z);
  free(nor_z2_x);
  free(nor_z2_z);

  return 0;
}

// strech grid base on arc length 
// use exponential function
int 
zt_arc_strech(gd_t *gdcurv, float coef)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t siz_icmp = gdcurv->siz_icmp;
  size_t iptr,iptr1,iptr2;
  float x_len,z_len,dh_len;
  float r, ratio, zeta;
  int n;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_temp; 
  float *z2d_temp;
  float *s;
  float *u;

  x2d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  z2d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");

 
  // line by line. i=0 -> i=nx-1 
  for(int i=0; i<nx; i++)
  {
    // copy old coords to temp space
    for(int k=0; k<nz; k++)
    {
      iptr1 = k*nx + i;     //(i,k)
      x2d_temp[k] = x2d[iptr1];
      z2d_temp[k] = z2d[iptr1];
    }
    // cal arc length
    for(int k=1; k<nz; k++)
    {
      x_len = x2d_temp[k] - x2d_temp[k-1];
      z_len = z2d_temp[k] - z2d_temp[k-1];
      dh_len = sqrt(pow(x_len,2) + pow(z_len,2));
      s[k] = s[k-1] + dh_len;
    }
    // arc length normalized
    for(int k=0; k<nz; k++)
    {
      u[k] = s[k]/s[nz-1];
    }
    for(int k=1; k<nz-1; k++)
    {
      zeta = (1.0*k)/(nz-1);
      r = single_exp(coef,zeta);
      for(int m=0; m<nz-1; m++)
      {
        if(r>=u[m] && r<u[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k*nx + i;
      x_len = x2d_temp[n+1] - x2d_temp[n];
      z_len = z2d_temp[n+1] - z2d_temp[n];
      ratio = (r - u[n])/(u[n+1]-u[n]);
      x2d[iptr] = x2d_temp[n] + x_len*ratio;
      z2d[iptr] = z2d_temp[n] + z_len*ratio;
    }
  }

  free(x2d_temp);
  free(z2d_temp);
  free(s);
  free(u);

  return 0;
}

// strech grid base on arc length 
// use exponential function
int 
xi_arc_strech(gd_t *gdcurv, float coef)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t siz_icmp = gdcurv->siz_icmp;
  size_t iptr,iptr1;
  float x_len,z_len,dh_len;
  float r, ratio, xi;
  int n;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_temp; 
  float *z2d_temp;
  float *s;
  float *u;

  x2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  z2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
 
  // line by line. k=0 -> k=nz-1 
  for(int k=0; k<nz; k++)
  {
    // copy old coords to temp space
    for(int i=0; i<nx; i++)
    {
      iptr1 = k*nx + i;     //(i,k)
      x2d_temp[i] = x2d[iptr1];
      z2d_temp[i] = z2d[iptr1];
    }
    // cal arc length
    for(int i=1; i<nx; i++)
    {
      x_len = x2d_temp[i] - x2d_temp[i-1];
      z_len = z2d_temp[i] - z2d_temp[i-1];
      dh_len = sqrt(pow(x_len,2) + pow(z_len,2));
      s[i] = s[i-1] + dh_len;
    }
    // arc length normalized
    for(int i=0; i<nx; i++)
    {
      u[i] = s[i]/s[nx-1];
    }
    for(int i=1; i<nx-1; i++)
    {
      xi = (1.0*i)/(nx-1);
      r = single_exp(coef,xi);
      for(int m=0; m<nx-1; m++)
      {
        if(r>=u[m] && r<u[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k*nx + i;
      x_len = x2d_temp[n+1] - x2d_temp[n];
      z_len = z2d_temp[n+1] - z2d_temp[n];
      ratio = (r - u[n])/(u[n+1]-u[n]);
      x2d[iptr] = x2d_temp[n] + x_len*ratio;
      z2d[iptr] = z2d_temp[n] + z_len*ratio;
    }
  }

  free(x2d_temp);
  free(z2d_temp);
  free(s);
  free(u);

  return 0;
}

// grid sample, linear interpolation  
int 
sample_interp(gd_t *gdcurv_new, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int nx_new = gdcurv_new->nx;
  int nz_new = gdcurv_new->nz;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_new = gdcurv_new->x2d;
  float *z2d_new = gdcurv_new->z2d;

  float *x2d_temp;
  float *z2d_temp;
  float *u;
  float *v;
  size_t iptr,iptr1,iptr2;
  float x_len,z_len;
  float r, ratio;
  int n;

  x2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  z2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nz, 0.0, "gd_curv_init");
  v = (float *)mem_calloc_1d_float(
              nx, 0.0, "gd_curv_init");

  // first interp zt direction 
  // line by line. i=0 -> i=nx-1 
  for(int i=0; i<nx; i++)
  {
    // point number normalized [0,1]
    for(int k=0; k<nz; k++)
    {
      u[k] = (1.0*k)/(nz-1);
    }

    for(int k_new=0; k_new<nz_new; k_new++)
    {
      r = (1.0*k_new)/(nz_new-1);
      for(int m=0; m<nz-1; m++)
      {
        if(r>=u[m] && r<u[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k_new*nx_new + i;
      iptr1 = n*nx + i;
      iptr2 = (n+1)*nx + i;
      x_len = x2d[iptr2] - x2d[iptr1];
      z_len = z2d[iptr2] - z2d[iptr1];
      ratio = (r - u[n])/(u[n+1]-u[n]);
      x2d_new[iptr] = x2d[iptr1] + x_len*ratio;
      z2d_new[iptr] = z2d[iptr1] + z_len*ratio;
    }
  }

  // then interp xi direction 
  // line by line. k=0 -> k=nz_new-1 
  for(int k_new=0; k_new<nz_new; k_new++)
  {
    // copy old coords to temp space
    for(int i=0; i<nx; i++)
    {
      iptr1 = k_new*nx_new + i;     //(i,k)
      x2d_temp[i] = x2d_new[iptr1];
      z2d_temp[i] = z2d_new[iptr1];
    }
    // point number normalized [0,1]
    for(int i=0; i<nx; i++)
    {
      v[i] = (1.0*i)/(nx-1);
    }

    for(int i_new=0; i_new<nx_new; i_new++)
    {
      r = (1.0*i_new)/(nx_new-1);
      for(int m=0; m<nx-1; m++)
      {
        if(r>=v[m] && r<v[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k_new*nx_new + i_new;
      x_len = x2d_temp[n+1] - x2d_temp[n];
      z_len = z2d_temp[n+1] - z2d_temp[n];
      ratio = (r - v[n])/(v[n+1]-v[n]);
      x2d_new[iptr] = x2d_temp[n] + x_len*ratio;
      z2d_new[iptr] = z2d_temp[n] + z_len*ratio;
    }
  }

  free(u);
  free(v);
  free(x2d_temp);
  free(z2d_temp);
  // use sample grid gdcurv_new, so free gdcurv space 
  free(gdcurv->v3d);

  return 0;
}

float 
single_exp(float coef, float zeta)
{
  float r;
  r = (exp(coef*zeta)-1)/(exp(coef)-1);
  return r; 
}
