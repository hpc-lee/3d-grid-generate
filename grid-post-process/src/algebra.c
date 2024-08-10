#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "algebra.h"
#include "lib_mem.h"

// stretch grid base on arc length 
int 
zt_arc_stretch(gd_t *gdcurv, float *arc_len)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t iptr,iptr1;
  float x_len,y_len,z_len,dh_len;
  float r, ratio, zt;
  int n;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *x3d_temp; 
  float *y3d_temp; 
  float *z3d_temp;
  float *s;
  float *u;

  x3d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  y3d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  z3d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");

 
  // line by line. 
  for(int j=0; j<ny; j++)
  {
    for(int i=0; i<nx; i++)
    {
      // copy old coords to temp space
      for(int k=0; k<nz; k++)
      {
        iptr1 = k*siz_iz + j*siz_iy + i;     //(i,j,k)
        x3d_temp[k] = x3d[iptr1];
        y3d_temp[k] = y3d[iptr1];
        z3d_temp[k] = z3d[iptr1];
      }
      // cal arc length
      for(int k=1; k<nz; k++)
      {
        x_len = x3d_temp[k] - x3d_temp[k-1];
        y_len = y3d_temp[k] - y3d_temp[k-1];
        z_len = z3d_temp[k] - z3d_temp[k-1];
        dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
        s[k] = s[k-1] + dh_len;
      }
      // arc length normalized
      for(int k=0; k<nz; k++)
      {
        u[k] = s[k]/s[nz-1];
      }
      for(int k=1; k<nz-1; k++)
      {
        zt = (1.0*k)/(nz-1);
        //r = single_exp(coef,zt);
        r = arc_len[k];
        for(int m=0; m<nz-1; m++)
        {
          if(r>=u[m] && r<u[m+1]) {
            n=m; 
            break;
          }
        }

        // linear interp
        iptr = k*siz_iz + j*siz_iy + i;
        x_len = x3d_temp[n+1] - x3d_temp[n];
        y_len = y3d_temp[n+1] - y3d_temp[n];
        z_len = z3d_temp[n+1] - z3d_temp[n];
        ratio = (r - u[n])/(u[n+1]-u[n]);
        x3d[iptr] = x3d_temp[n] + x_len*ratio;
        y3d[iptr] = y3d_temp[n] + y_len*ratio;
        z3d[iptr] = z3d_temp[n] + z_len*ratio;
      }
    }
  }

  free(x3d_temp);
  free(y3d_temp);
  free(z3d_temp);
  free(s);
  free(u);

  return 0;
}

// stretch grid base on arc length 
int 
et_arc_stretch(gd_t *gdcurv, float *arc_len)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t iptr,iptr1;
  float x_len,y_len,z_len,dh_len;
  float r, ratio, et;
  int n;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *x3d_temp; 
  float *y3d_temp; 
  float *z3d_temp;
  float *s;
  float *u;

  x3d_temp = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  y3d_temp = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  z3d_temp = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
 
  // line by line.
  for(int k=0; k<nz; k++)
  {
    for(int i=0; i<nx; i++)
    {
      // copy old coords to temp space
      for(int j=0; j<ny; j++)
      {
        iptr1 = k*siz_iz + j*siz_iy + i;     //(i,j,k)
        x3d_temp[j] = x3d[iptr1];
        y3d_temp[j] = y3d[iptr1];
        z3d_temp[j] = z3d[iptr1];
      }
      // cal arc length
      for(int j=1; j<ny; j++)
      {
        x_len = x3d_temp[j] - x3d_temp[j-1];
        y_len = y3d_temp[j] - y3d_temp[j-1];
        z_len = z3d_temp[j] - z3d_temp[j-1];
        dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
        s[j] = s[j-1] + dh_len;
      }
      // arc length normalized
      for(int j=0; j<ny; j++)
      {
        u[j] = s[j]/s[ny-1];
      }
      for(int j=1; j<ny-1; j++)
      {
        et = (1.0*j)/(ny-1);
        //r = single_exp(coef,et);
        r = arc_len[j];
        for(int m=0; m<ny-1; m++)
        {
          if(r>=u[m] && r<u[m+1]) {
            n=m; 
            break;
          }
        }

        // linear interp
        iptr = k*siz_iz + j*siz_iy + i;
        x_len = x3d_temp[n+1] - x3d_temp[n];
        y_len = y3d_temp[n+1] - y3d_temp[n];
        z_len = z3d_temp[n+1] - z3d_temp[n];
        ratio = (r - u[n])/(u[n+1]-u[n]);
        x3d[iptr] = x3d_temp[n] + x_len*ratio;
        y3d[iptr] = y3d_temp[n] + y_len*ratio;
        z3d[iptr] = z3d_temp[n] + z_len*ratio;
      }
    }
  }

  free(x3d_temp);
  free(y3d_temp);
  free(z3d_temp);
  free(s);
  free(u);

  return 0;
}

// stretch grid base on arc length 
int 
xi_arc_stretch(gd_t *gdcurv, float *arc_len)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t iptr,iptr1;
  float x_len,y_len,z_len,dh_len;
  float r, ratio, xi;
  int n;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *x3d_temp; 
  float *y3d_temp; 
  float *z3d_temp;
  float *s;
  float *u;

  x3d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  y3d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  z3d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
 
  // line by line.
  for(int k=0; k<nz; k++)
  {
    for(int j=0; j<ny; j++)
    {
      // copy old coords to temp space
      for(int i=0; i<nx; i++)
      {
        iptr1 = k*siz_iz + j*siz_iy + i;     //(i,j,k)
        x3d_temp[i] = x3d[iptr1];
        y3d_temp[i] = y3d[iptr1];
        z3d_temp[i] = z3d[iptr1];
      }
      // cal arc length
      for(int i=1; i<nx; i++)
      {
        x_len = x3d_temp[i] - x3d_temp[i-1];
        y_len = y3d_temp[i] - y3d_temp[i-1];
        z_len = z3d_temp[i] - z3d_temp[i-1];
        dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
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
        //r = single_exp(coef,xi);
        r = arc_len[i];
        for(int m=0; m<nx-1; m++)
        {
          if(r>=u[m] && r<u[m+1]) {
            n=m; 
            break;
          }
        }

        // linear interp
        iptr = k*siz_iz + j*siz_iy + i;
        x_len = x3d_temp[n+1] - x3d_temp[n];
        y_len = y3d_temp[n+1] - y3d_temp[n];
        z_len = z3d_temp[n+1] - z3d_temp[n];
        ratio = (r - u[n])/(u[n+1]-u[n]);
        x3d[iptr] = x3d_temp[n] + x_len*ratio;
        y3d[iptr] = y3d_temp[n] + y_len*ratio;
        z3d[iptr] = z3d_temp[n] + z_len*ratio;
      }
    }
  }

  free(x3d_temp);
  free(y3d_temp);
  free(z3d_temp);
  free(s);
  free(u);

  return 0;
}

// grid sample, linear interpolation  
int 
sample_interp(gd_t *gdcurv_new, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int nx_new = gdcurv_new->nx;
  int ny_new = gdcurv_new->ny;
  int nz_new = gdcurv_new->nz;
  size_t siz_iy_new = gdcurv_new->siz_iy;
  size_t siz_iz_new = gdcurv_new->siz_iz;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *x3d_new = gdcurv_new->x3d;
  float *y3d_new = gdcurv_new->y3d;
  float *z3d_new = gdcurv_new->z3d;

  float *x3d_temp_x, *x3d_temp_y;
  float *y3d_temp_x, *y3d_temp_y;
  float *z3d_temp_x, *z3d_temp_y;
  float *u_x, *u_y, *u_z;
  size_t iptr,iptr1,iptr2;
  float x_len,y_len,z_len;
  float r, ratio;
  int n;

  x3d_temp_x = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  y3d_temp_x = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  z3d_temp_x = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  x3d_temp_y = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  y3d_temp_y = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  z3d_temp_y = (float *)mem_calloc_1d_float(
              ny, 0.0, "init");
  u_z = (float *)mem_calloc_1d_float(
              nz, 0.0, "gd_curv_init");
  u_y = (float *)mem_calloc_1d_float(
              ny, 0.0, "gd_curv_init");
  u_x = (float *)mem_calloc_1d_float(
              nx, 0.0, "gd_curv_init");

  // first interp zt direction 
  for(int j=0; j<ny; j++)
  {
    for(int i=0; i<nx; i++)
    {
      // point number normalized [0,1]
      for(int k=0; k<nz; k++)
      {
        u_z[k] = (1.0*k)/(nz-1);
      }

      for(int k1=0; k1<nz_new; k1++)
      {
        r = (1.0*k1)/(nz_new-1);
        for(int m=0; m<nz-1; m++)
        {
          if(r>=u_z[m] && r<u_z[m+1]) {
            n=m; 
            break;
          }
        }

        // linear interp
        iptr = k1*siz_iz_new + j*siz_iy_new + i;
        iptr1 = n*siz_iz + j*siz_iy + i;
        iptr2 = (n+1)*siz_iz + j*siz_iy + i;
        x_len = x3d[iptr2] - x3d[iptr1];
        y_len = y3d[iptr2] - y3d[iptr1];
        z_len = z3d[iptr2] - z3d[iptr1];
        ratio = (r - u_z[n])/(u_z[n+1]-u_z[n]);
        x3d_new[iptr] = x3d[iptr1] + x_len*ratio;
        y3d_new[iptr] = y3d[iptr1] + y_len*ratio;
        z3d_new[iptr] = z3d[iptr1] + z_len*ratio;
      }
    }
  }

  // then interp et direction 
  for(int k1=0; k1<nz_new; k1++)
  {
    for(int i=0; i<nx; i++)
    {
      // copy old coords to temp space
      for(int j=0; j<ny; j++)
      {
        iptr = k1*siz_iz_new + j*siz_iy_new + i;
        x3d_temp_y[j] = x3d_new[iptr];
        y3d_temp_y[j] = y3d_new[iptr];
        z3d_temp_y[j] = z3d_new[iptr];
      }
      // point number normalized [0,1]
      for(int j=0; j<ny; j++)
      {
        u_y[j] = (1.0*j)/(ny-1);
      }

      for(int j1=0; j1<ny_new; j1++)
      {
        r = (1.0*j1)/(ny_new-1);
        for(int m=0; m<ny-1; m++)
        {
          if(r>=u_y[m] && r<u_y[m+1]) {
            n=m; 
            break;
          }
        }

        // linear interp
        iptr = k1*siz_iz_new + j1*siz_iy_new + i;
        x_len = x3d_temp_y[n+1] - x3d_temp_y[n];
        y_len = y3d_temp_y[n+1] - y3d_temp_y[n];
        z_len = z3d_temp_y[n+1] - z3d_temp_y[n];
        ratio = (r - u_y[n])/(u_y[n+1]-u_y[n]);
        x3d_new[iptr] = x3d_temp_y[n] + x_len*ratio;
        y3d_new[iptr] = y3d_temp_y[n] + y_len*ratio;
        z3d_new[iptr] = z3d_temp_y[n] + z_len*ratio;
      }
    }
  }

  // finally interp xi direction 
  for(int k1=0; k1<nz_new; k1++)
  {
    for(int j1=0; j1<ny_new; j1++)
    {
      for(int i=0; i<nx; i++)
      // copy old coords to temp space
      {
        iptr = k1*siz_iz_new + j1*siz_iy_new + i;
        x3d_temp_x[i] = x3d_new[iptr];
        y3d_temp_x[i] = y3d_new[iptr];
        z3d_temp_x[i] = z3d_new[iptr];
      }
      // point number normalized [0,1]
      for(int i=0; i<nx; i++)
      {
        u_x[i] = (1.0*i)/(nx-1);
      }

      for(int i1=0; i1<nx_new; i1++)
      {
        r = (1.0*i1)/(nx_new-1);
        for(int m=0; m<nx-1; m++)
        {
          if(r>=u_x[m] && r<u_x[m+1]) {
            n=m; 
            break;
          }
        }

        // linear interp
        iptr = k1*siz_iz_new + j1*siz_iy_new + i1;
        x_len = x3d_temp_x[n+1] - x3d_temp_x[n];
        y_len = y3d_temp_x[n+1] - y3d_temp_x[n];
        z_len = z3d_temp_x[n+1] - z3d_temp_x[n];
        ratio = (r - u_x[n])/(u_x[n+1]-u_x[n]);
        x3d_new[iptr] = x3d_temp_x[n] + x_len*ratio;
        y3d_new[iptr] = y3d_temp_x[n] + y_len*ratio;
        z3d_new[iptr] = z3d_temp_x[n] + z_len*ratio;
      }
    }
  }

  free(u_x);
  free(u_y);
  free(u_z);
  free(x3d_temp_x);
  free(y3d_temp_x);
  free(z3d_temp_x);
  free(x3d_temp_y);
  free(y3d_temp_y);
  free(z3d_temp_y);
  // use sample grid gdcurv_new, so free gdcurv space 
  free(gdcurv->v4d);

  return 0;
}

float 
single_exp(float coef, float zeta)
{
  float r;
  r = (exp(coef*zeta)-1)/(exp(coef)-1);
  return r; 
}
