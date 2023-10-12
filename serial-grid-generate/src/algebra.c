#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "algebra.h"
#include "lib_mem.h"

// linear tfi interpolation
// U means x-direction
// V means y-direction
// W means z-direction
int linear_tfi(gd_t *gdcurv)
{
  int nx = gdcurv->nx; 
  int ny = gdcurv->ny; 
  int nz = gdcurv->nz; 
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float xi,et,zt;
  float a0,a1,b0,b1,c0,c1;
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  size_t iptr5,iptr6,iptr7,iptr8;
  size_t iptr9,iptr10,iptr11,iptr12;
  float U_x,V_x,W_x,UV_x,UW_x,VW_x,UVW_x;
  float U_y,V_y,W_y,UV_y,UW_y,VW_y,UVW_y;
  float U_z,V_z,W_z,UV_z,UW_z,VW_z,UVW_z;
  
  for (int k=1; k<nz-1; k++) {
    for (int j=1; j<ny-1; j++) {
      for (int i=1; i<nx-1; i++)
      {
        xi = (1.0*i)/(nx-1); // 1.0* equal int to float
        et = (1.0*j)/(ny-1);
        zt = (1.0*k)/(nz-1);
        a0 = 1-xi;
        a1 = xi;
        b0 = 1-et;
        b1 = et;
        c0 = 1-zt;
        c1 = zt;

        iptr = k*siz_iz + j*siz_iy + i;  //(k,j,i)
        // 6 face points
        iptr1 = k*siz_iz + j*siz_iy;     //(k,j,0)
        iptr2 = k*siz_iz + j*siz_iy + (nx-1); //(k,j,nx-1)
        iptr3 = k*siz_iz + i;         //(k,0,i)
        iptr4 = k*siz_iz + (ny-1)*siz_iy + i; //(k,ny-1,i)
        iptr5 = j*siz_iy + i; //(0,j,i)
        iptr6 = (nz-1)*siz_iz + j*siz_iy + i; //(nz-1,j,i)
        U_x = a0*x3d[iptr1] + a1*x3d[iptr2];
        V_x = b0*x3d[iptr3] + b1*x3d[iptr4];
        W_x = c0*x3d[iptr5] + c1*x3d[iptr6];

        U_y = a0*y3d[iptr1] + a1*y3d[iptr2];
        V_y = b0*y3d[iptr3] + b1*y3d[iptr4];
        W_y = c0*y3d[iptr5] + c1*y3d[iptr6];

        U_z = a0*z3d[iptr1] + a1*z3d[iptr2];
        V_z = b0*z3d[iptr3] + b1*z3d[iptr4];
        W_z = c0*z3d[iptr5] + c1*z3d[iptr6];
        // 12 edge points
        iptr1 = k*siz_iz;   //(k,0,0)
        iptr2 = k*siz_iz + (ny-1)*siz_iy;   //(k,ny-1,0)
        iptr3 = k*siz_iz + nx-1;   //(k,0,nx-1)
        iptr4 = k*siz_iz + (ny-1)*siz_iy + (nx-1);   //(k,ny-1,nx-1)
        iptr5 = j*siz_iy;  //(0,j,0)
        iptr6 = (nz-1)*siz_iz + j*siz_iy;  //(nz-1,j,0)
        iptr7 = j*siz_iy + (nx-1);  //(0,j,nx-1)
        iptr8 = (nz-1)*siz_iz + j*siz_iy + (nx-1);  //(nz-1,j,nx)
        iptr9 = i; //(0,0,i)
        iptr10 = (nz-1)*siz_iz + i; //(nz-1,0,i)
        iptr11 = (ny-1)*siz_iy + i; //(0,ny-1,i)
        iptr12 = (nz-1)*siz_iz + (ny-1)*siz_iy + i; //(nz-1,ny-1,i)

        UV_x = a0*b0*x3d[iptr1]  + a0*b1*x3d[iptr2]
             + a1*b0*x3d[iptr3]  + a1*b1*x3d[iptr4];
        UW_x = a0*c0*x3d[iptr5]  + a0*c1*x3d[iptr6]
             + a1*c0*x3d[iptr7]  + a1*c1*x3d[iptr8];
        VW_x = b0*c0*x3d[iptr9]  + b0*c1*x3d[iptr10]
             + b1*c0*x3d[iptr11] + b1*c1*x3d[iptr12];

        UV_y = a0*b0*y3d[iptr1]  + a0*b1*y3d[iptr2]
             + a1*b0*y3d[iptr3]  + a1*b1*y3d[iptr4];
        UW_y = a0*c0*y3d[iptr5]  + a0*c1*y3d[iptr6]
             + a1*c0*y3d[iptr7]  + a1*c1*y3d[iptr8];
        VW_y = b0*c0*y3d[iptr9]  + b0*c1*y3d[iptr10]
             + b1*c0*y3d[iptr11] + b1*c1*y3d[iptr12];

        UV_z = a0*b0*z3d[iptr1]  + a0*b1*z3d[iptr2]
             + a1*b0*z3d[iptr3]  + a1*b1*z3d[iptr4];
        UW_z = a0*c0*z3d[iptr5]  + a0*c1*z3d[iptr6]
             + a1*c0*z3d[iptr7]  + a1*c1*z3d[iptr8];
        VW_z = b0*c0*z3d[iptr9]  + b0*c1*z3d[iptr10]
             + b1*c0*z3d[iptr11] + b1*c1*z3d[iptr12];

        // 8 corner point
        iptr1 = 0; //(0,0,0)
        iptr2 = nx-1; //(0,0,nx-1)
        iptr3 = (ny-1)*siz_iy + (nx-1); //(0,ny-1,nx-1)
        iptr4 = (ny-1)*siz_iy; //(0,ny-1,0)
        iptr5 = (nz-1)*siz_iz; //(nz-1,0,0)
        iptr6 = (nz-1)*siz_iz + (nx-1); //(nz-1,0,nx-1)
        iptr7 = (nz-1)*siz_iz + (ny-1)*siz_iy + (nx-1); //(nz-1,ny-1,nx-1)
        iptr8 = (nz-1)*siz_iz + (ny-1)*siz_iy; //(nz-1,ny-1,0)
        UVW_x = a0*b0*c0*x3d[iptr1] + a1*b0*c0*x3d[iptr2]
               +a1*b1*c0*x3d[iptr3] + a0*b1*c0*x3d[iptr4]
               +a0*b0*c1*x3d[iptr5] + a1*b0*c1*x3d[iptr6]
               +a1*b1*c1*x3d[iptr7] + a0*b1*c1*x3d[iptr8];

        UVW_y = a0*b0*c0*y3d[iptr1] + a1*b0*c0*y3d[iptr2]
               +a1*b1*c0*y3d[iptr3] + a0*b1*c0*y3d[iptr4]
               +a0*b0*c1*y3d[iptr5] + a1*b0*c1*y3d[iptr6]
               +a1*b1*c1*y3d[iptr7] + a0*b1*c1*y3d[iptr8];

        UVW_z = a0*b0*c0*z3d[iptr1] + a1*b0*c0*z3d[iptr2]
               +a1*b1*c0*z3d[iptr3] + a0*b1*c0*z3d[iptr4]
               +a0*b0*c1*z3d[iptr5] + a1*b0*c1*z3d[iptr6]
               +a1*b1*c1*z3d[iptr7] + a0*b1*c1*z3d[iptr8];

        x3d[iptr] = U_x + V_x + W_x - UV_x - UW_x - VW_x + UVW_x;
        y3d[iptr] = U_y + V_y + W_y - UV_y - UW_y - VW_y + UVW_y;
        z3d[iptr] = U_z + V_z + W_z - UV_z - UW_z - VW_z + UVW_z;
      }
    }
  }
  
  return 0;
}

