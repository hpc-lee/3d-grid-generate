#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"
#include "lib_math.h"
#include "lib_mem.h"
#include "constants.h"

// solve tridiagonal linear system equaltion
// using thomas method
// Due to the sparse coefficient matrix of the tridiagonal equation,
// the computational complexity is proportional to n,
// rather than the n^3 of Gaussian elimination
//
// [b1 c1        ]
// |a2 b2 c2     |
// |   a3 b3 c3  |
// [             ]
//
// a = element of lower diagonal
// b = element of main diagonal
// c = element of up diagonal
// d = right hand item
// u is unknow vector
int
thomas(int n, float *a, float *b, float *c, float *d_x, 
       float *d_y, float *d_z, float *u_x, float *u_y, float *u_z)
{
  float factor;
  for(int i=1; i<n; i++)
  {
    factor = a[i]/b[i-1];
    b[i] = b[i] - factor*c[i-1];
    d_x[i] = d_x[i] - factor*d_x[i-1];
    d_y[i] = d_y[i] - factor*d_y[i-1];
    d_z[i] = d_z[i] - factor*d_z[i-1];
  }

  u_x[n-1] = d_x[n-1]/b[n-1];
  u_y[n-1] = d_y[n-1]/b[n-1];
  u_z[n-1] = d_z[n-1]/b[n-1];
  for(int i=n-2; i>=0; i--)
  {
    u_x[i] = (d_x[i]-c[i]*u_x[i+1])/b[i];
    u_y[i] = (d_y[i]-c[i]*u_y[i+1])/b[i];
    u_z[i] = (d_z[i]-c[i]*u_z[i+1])/b[i];
  }

  return 0;
}

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

/*
  solve block tridiagonal linear system equaltion
  using thomas method
  Due to the sparse coefficient matrix of the tridiagonal equation,
  the computational complexity is proportional to n,
  rather than the n^3 of Gaussian elimination 

  [b1 c1        ]
  |a2 b2 c2     |
  |   a3 b3 c3  |
  [             ]

  a = square matrix element of lower diagonal
  b = square matrix element of main diagonal
  c = square matrix element of up diagonal
  d = right hand item, each element size is n*1

    [G1 0  0       ]
    |a2 G2 0       |
    |   a3 G3      |
  L=|              |
    |              |
    |              |
    [         an Gn]
                    
    [I1 D1 0            ]
    |0  I2 D2           |
    |   0  I3           |
  U=|                   |
    |                   |
    |          In-1 Dn-1|
    [          0    In  ]

  LU*x = d
  Ly=d
  Ux=y
*/
int
thomas_block(int n, double *a, double *b, double *c, double *d,
             double *xz, double *D, double *y)
{
  double mat_a[2][2], mat_b[2][2], mat_c[2][2], vec_d[2];
  double mat_G[2][2], mat_D[2][2], vec_y[2], vec_xz[2];
  double det, mat1[2][2], vec1[2], vec2[2];
  size_t iptr1,iptr2,iptr3,iptr4;

  // i=0
  for(int ii=0; ii<2; ii++) {
    for(int jj=0; jj<2; jj++) {
      mat_G[ii][jj] = b[ii*2+jj];
      mat_c[ii][jj] = c[ii*2+jj];
    }
    vec_d[ii] = d[ii];
  }

  det = mat_G[0][0]*mat_G[1][1] - mat_G[0][1]*mat_G[1][0];
  if(fabs(det)<0.000001)
  {
    fprintf(stderr,"b{1} is singular, stop calculate");
    exit(1);
  } 
  mat_invert2x2(mat_G);
  mat_mul2x2(mat_G,mat_c,mat_D);
  mat_mul2x1(mat_G,vec_d,vec_y);
  for(int ii=0; ii<2; ii++) { 
    for(int jj=0; jj<2; jj++) {
      D[ii*2+jj] = mat_D[ii][jj];
    }
    y[ii] = vec_y[ii];
  }

  // i=1~n-1
  for(int i=1; i<n; i++)
  {
    iptr1 = i*CONST_NDIM*CONST_NDIM;
    iptr2 = i*CONST_NDIM;
    iptr3 = (i-1)*CONST_NDIM*CONST_NDIM;
    iptr4 = (i-1)*CONST_NDIM;
    for(int ii=0; ii<2; ii++) { 
      for(int jj=0; jj<2; jj++) {
        mat_a[ii][jj] = a[iptr1+ii*2+jj];
        mat_b[ii][jj] = b[iptr1+ii*2+jj];
        mat_c[ii][jj] = c[iptr1+ii*2+jj];
        mat_D[ii][jj] = D[iptr3+ii*2+jj];
      }
      vec_d[ii] = d[iptr2+ii];
      vec_y[ii] = y[iptr4+ii];
    }
    mat_mul2x2(mat_a,mat_D,mat1);  // a(i)*D(i-1)
    mat_sub2x2(mat_b,mat1,mat_G);  // G(i) = b(i)-a(i)*D(i-1)
    mat_invert2x2(mat_G);          // inv(G(i))
    mat_mul2x2(mat_G,mat_c,mat_D); // D(i) = inv(G(i))*c(i)
    mat_mul2x1(mat_a,vec_y,vec1);  // a(i) * y(i-1)
    vec_sub2x1(vec_d,vec1,vec2);   // d(i) - a(i) * y(i-1)
    mat_mul2x1(mat_G,vec2,vec_y);  // y(i) = inv(G(i))*(d(i) - a(i) * y(i-1))

    for(int ii=0; ii<2; ii++) { 
      for(int jj=0; jj<2; jj++) {
        D[iptr1+ii*2+jj] = mat_D[ii][jj]; 
      }
      y[iptr2+ii] = vec_y[ii];
    }
  }
  
  // i=n-1
  iptr2 = (n-1)*CONST_NDIM;
  for(int ii=0; ii<2; ii++) { 
    xz[iptr2+ii] = y[iptr2+ii];
  }

  for(int i=n-2; i>=0; i--)
  {
    iptr1 = i*CONST_NDIM*CONST_NDIM;
    iptr2 = i*CONST_NDIM;
    iptr4 = (i+1)*CONST_NDIM;
    
    for(int ii=0; ii<2; ii++) { 
      for(int jj=0; jj<2; jj++) {
        mat_D[ii][jj] = D[iptr1+ii*2+jj]; 
      }
      vec_xz[ii] = xz[iptr4+ii];
    }
    mat_mul2x1(mat_D,vec_xz,vec1);  // D(i)*xz(i+1)
    for(int ii=0; ii<2; ii++) 
    { 
      xz[iptr2+ii] = y[iptr2+ii] - vec1[ii]; // xz(i) = y(i) - D(i)*xz(i+1) 
    }
  }
  
  return 0;
}

