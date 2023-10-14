#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"
#include "lib_math.h"
#include "lib_mem.h"
#include "constants.h"

int
update_SOR(float *x3d, float *y3d, float *z3d, float *x3d_tmp,
           float *y3d_tmp, float *z3d_tmp, int nx, int ny, int nz,
           float *P, float *Q, float *R, float w)
{
  size_t iptr,iptr1,iptr2,iptr3,iptr4,iptr5,iptr6;
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xiet,y_xiet,z_xiet;
  float x_xizt,y_xizt,z_xizt;
  float x_etzt,y_etzt,z_etzt;
  float g11,g22,g33,g12,g13,g23,coef;
  float alpha1,alpha2,alpha3;
  float beta12,beta23,beta13;

  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = k*siz_iz + j*siz_iy + (i+1);  // (i+1,j,k)
        iptr2 = k*siz_iz + j*siz_iy + (i-1);  // (i-1,j,k)
        iptr3 = k*siz_iz + (j+1)*siz_iy + i;  // (i,j+1,k)
        iptr4 = k*siz_iz + (j-1)*siz_iy + i;  // (i,j-1,k)
        iptr5 = (k+1)*siz_iz + j*siz_iy + i;  // (i,j,k+1)
        iptr6 = (k-1)*siz_iz + j*siz_iy + i;  // (i,j,k-1)

        x_xi = 0.5*(x3d[iptr1] - x3d_tmp[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d_tmp[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d_tmp[iptr2]);

        x_et = 0.5*(x3d[iptr3] - x3d_tmp[iptr4]);
        y_et = 0.5*(y3d[iptr3] - y3d_tmp[iptr4]);
        z_et = 0.5*(z3d[iptr3] - z3d_tmp[iptr4]);

        x_zt = 0.5*(x3d[iptr5] - x3d_tmp[iptr6]);
        y_zt = 0.5*(y3d[iptr5] - y3d_tmp[iptr6]);
        z_zt = 0.5*(z3d[iptr5] - z3d_tmp[iptr6]);

        iptr1 = k*siz_iz + (j-1)*siz_iy + (i-1);   // (i-1,j-1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + (i+1);   // (i+1,j-1,k)
        iptr3 = k*siz_iz + (j+1)*siz_iy + (i+1);   // (i+1,j+1,k)
        iptr4 = k*siz_iz + (j+1)*siz_iy + (i-1);   // (i-1,j+1,k)

        x_xiet = 0.25*(x3d[iptr3] + x3d_tmp[iptr1] - x3d_tmp[iptr2] - x3d[iptr4]);
        y_xiet = 0.25*(y3d[iptr3] + y3d_tmp[iptr1] - y3d_tmp[iptr2] - y3d[iptr4]);
        z_xiet = 0.25*(z3d[iptr3] + z3d_tmp[iptr1] - z3d_tmp[iptr2] - z3d[iptr4]);

        iptr1 = (k-1)*siz_iz + j*siz_iy + (i-1);   // (i-1,j,k-1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + (i+1);   // (i+1,j,k-1)
        iptr3 = (k+1)*siz_iz + j*siz_iy + (i+1);   // (i+1,j,k+1)
        iptr4 = (k+1)*siz_iz + j*siz_iy + (i-1);   // (i-1,j,k+1)

        x_xizt = 0.25*(x3d[iptr3] + x3d_tmp[iptr1] - x3d_tmp[iptr2] - x3d[iptr4]);
        y_xizt = 0.25*(y3d[iptr3] + y3d_tmp[iptr1] - y3d_tmp[iptr2] - y3d[iptr4]);
        z_xizt = 0.25*(z3d[iptr3] + z3d_tmp[iptr1] - z3d_tmp[iptr2] - z3d[iptr4]);

        iptr1 = (k-1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-1)
        iptr2 = (k-1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-1)
        iptr3 = (k+1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k+1)
        iptr4 = (k+1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k+1)

        x_etzt = 0.25*(x3d[iptr3] + x3d_tmp[iptr1] - x3d_tmp[iptr2] - x3d[iptr4]);
        y_etzt = 0.25*(y3d[iptr3] + y3d_tmp[iptr1] - y3d_tmp[iptr2] - y3d[iptr4]);
        z_etzt = 0.25*(z3d[iptr3] + z3d_tmp[iptr1] - z3d_tmp[iptr2] - z3d[iptr4]);


        g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
        g22 = x_et*x_et + y_et*y_et + z_et*z_et;
        g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
        g12 = x_xi*x_et + y_xi*y_et + z_xi*z_et;
        g13 = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
        g23 = x_et*x_zt + y_et*y_zt + z_et*z_zt;

        alpha1 = g22*g33 - g23*g23;
        alpha2 = g11*g33 - g13*g13;
        alpha3 = g11*g22 - g12*g12;
        beta12 = g13*g23 - g12*g33;
        beta23 = g12*g13 - g11*g23;
        beta13 = g12*g23 - g13*g22;

        coef = 0.5/(alpha1 + alpha2 + alpha3);
        
        iptr  = k*siz_iz + j*siz_iy + i;     // (i,j,k)
        iptr1 = k*siz_iz + j*siz_iy + (i+1); // (i+1,j,k)
        iptr2 = k*siz_iz + j*siz_iy + (i-1); // (i-1,j,k)
        iptr3 = k*siz_iz + (j+1)*siz_iy + i;  // (i,j+1,k)
        iptr4 = k*siz_iz + (j-1)*siz_iy + i;  // (i,j-1,k)
        iptr5 = (k+1)*siz_iz + j*siz_iy + i;  // (i,j,k+1)
        iptr6 = (k-1)*siz_iz + j*siz_iy + i;  // (i,j,k-1)

        x3d_tmp[iptr] = coef*(alpha1*(x3d[iptr1]+x3d_tmp[iptr2]) + alpha2*(x3d[iptr3]
                        + x3d_tmp[iptr4]) + alpha3*(x3d[iptr5]+x3d_tmp[iptr6])
                        + 2*(beta12*x_xiet + beta23*x_etzt + beta13*x_xizt)
                        + alpha1*P[iptr]*x_xi + alpha2*Q[iptr]*x_et + alpha3*R[iptr]*x_zt);

        x3d_tmp[iptr] = w*x3d_tmp[iptr] + (1-w)*x3d[iptr];

        y3d_tmp[iptr] = coef*(alpha1*(y3d[iptr1]+y3d_tmp[iptr2]) + alpha2*(y3d[iptr3]
                        + y3d_tmp[iptr4]) + alpha3*(y3d[iptr5]+y3d_tmp[iptr6])
                        + 2*(beta12*y_xiet + beta23*y_etzt + beta13*y_xizt)
                        + alpha1*P[iptr]*y_xi + alpha2*Q[iptr]*y_et + alpha3*R[iptr]*y_zt);

        y3d_tmp[iptr] = w*y3d_tmp[iptr] + (1-w)*y3d[iptr];

        z3d_tmp[iptr] = coef*(alpha1*(z3d[iptr1]+z3d_tmp[iptr2]) + alpha2*(z3d[iptr3]
                        + z3d_tmp[iptr4]) + alpha3*(z3d[iptr5]+z3d_tmp[iptr6])
                        + 2*(beta12*z_xiet + beta23*z_etzt + beta13*z_xizt)
                        + alpha1*P[iptr]*z_xi + alpha2*Q[iptr]*z_et + alpha3*R[iptr]*z_zt);

        z3d_tmp[iptr] = w*z3d_tmp[iptr] + (1-w)*z3d[iptr];

      }
    }
  }

  return 0;
}
