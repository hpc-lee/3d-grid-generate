#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "hyperbolic.h"
#include "solver.h"
#include "lib_mem.h"
#include "lib_math.h"
#include "constants.h"

int 
hyper_gene(gd_t *gdcurv, float coef, int o2i, int bdry_itype, float epsilon)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int n = nx-2;  // not include bdry 2 points

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *step = gdcurv->step;

  float *coef_e = NULL;
  coef_e = (float *)mem_calloc_1d_float(n, 0.0, "init");
  float *area = NULL;
  area = (float *)mem_calloc_1d_float(nx*2, 0.0, "init");
  // malloc space for thomas_block method 
  double *a = NULL;
  a = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *b = NULL;
  b = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *c = NULL;
  c = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *d = NULL;
  d = (double *)mem_calloc_1d_double(n*CONST_NDIM, 0.0, "init");
  double *xz = NULL;
  xz = (double *)mem_calloc_1d_double(n*CONST_NDIM, 0.0, "init");
  double *D = NULL;
  D = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *y = NULL;
  y = (double *)mem_calloc_1d_double(n*CONST_NDIM, 0.0, "init");

  for(int k=1; k<nz; k++)
  {
    cal_smooth_coef(coef,x2d,z2d,nx,nz,k,coef_e);
    cal_matrix(x2d,z2d,nx,k,step,a,b,c,d,area);
    modify_smooth(x2d,z2d,nx,k,a,b,c,d,coef_e);
    modify_bdry(n,a,b,c,d,epsilon,bdry_itype);
    thomas_block(n,a,b,c,d,xz,D,y);
    assign_coords(xz,x2d,z2d,nx,k,epsilon,bdry_itype);
  }

  if(o2i == 1)
  {
    fprintf(stdout,"hyperbolic method, inner bdry(k=0), outer bdry(nz-1)\n");
    fprintf(stdout,"we default set read init bdry is inner bdry(k=0)\n");
    fprintf(stdout,"so if the init bdry is outer bdry actually, must be flip\n");
    flip_coord(x2d,nx,nz);
    flip_coord(z2d,nx,nz);
  }

  free(coef_e);
  free(area);
  free(a);
  free(b);
  free(c);
  free(d);
  free(xz);
  free(D);
  free(y);

  return 0;
}

int
cal_smooth_coef(float coef, float *x2d, float *z2d, int nx, int nz, int k, float *coef_e)
{
  float S;
  size_t iptr1, iptr2, iptr3, iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float xi_len,zt_len,N_xi;
  float x_xi_plus,z_xi_plus,x_xi_minus,z_xi_minus;
  float xi_plus1,xi_minus1,xi_plus2,xi_minus2;
  float d1,d2,delta,delta_mdfy;
  float x_plus,z_plus,x_minus,z_minus;
  float dot,det,theta,alpha;

  if(nz<=20) {
    S = sqrt((1.0*k)/(nz-1));
  }

  if(nz>20) 
  {
    if(k<20)
    {
      S = sqrt((1.0*k)/19);
    } else {
      S = 1.0;
    }
  }

  // k=1 do nothing
  if(k>1)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i+1;   // (i+1,k-1)
      iptr2 = (k-1)*nx + i-1;   // (i-1,k-1)
      iptr3 = (k-1)*nx + i;     // (i,k-1)
      iptr4 = (k-2)*nx + i;     // (i,k-2)
      x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]); 
      x_zt = x2d[iptr3] - x2d[iptr4];
      z_zt = z2d[iptr3] - z2d[iptr4];
      xi_len = sqrt(pow(x_xi,2) + pow(z_xi,2));
      zt_len = sqrt(pow(x_zt,2) + pow(z_zt,2));
      N_xi = zt_len/xi_len;

      iptr1 = (k-2)*nx + i+1;   // (i+1,k-2)
      iptr2 = (k-2)*nx + i;     // (i,  k-2)
      iptr3 = (k-2)*nx + i-1;   // (i-1,k-2)
      x_xi_plus = x2d[iptr1] - x2d[iptr2];
      z_xi_plus = z2d[iptr1] - z2d[iptr2];
      xi_plus1 = sqrt(pow(x_xi_plus,2) + pow(z_xi_plus,2));
      x_xi_minus = x2d[iptr3] - x2d[iptr2];
      z_xi_minus = z2d[iptr3] - z2d[iptr2];
      xi_minus1 = sqrt(pow(x_xi_minus,2) + pow(z_xi_minus,2));

      iptr1 = (k-1)*nx + i+1;   // (i+1,k-1)
      iptr2 = (k-1)*nx + i;     // (i,  k-1)
      iptr3 = (k-1)*nx + i-1;   // (i-1,k-1)
      x_xi_plus = x2d[iptr1] - x2d[iptr2];
      z_xi_plus = z2d[iptr1] - z2d[iptr2];
      xi_plus2 = sqrt(pow(x_xi_plus,2) + pow(z_xi_plus,2));
      x_xi_minus = x2d[iptr3] - x2d[iptr2];
      z_xi_minus = z2d[iptr3] - z2d[iptr2];
      xi_minus2 = sqrt(pow(x_xi_minus,2) + pow(z_xi_minus,2));

      d1 = xi_plus1 + xi_minus1;
      d2 = xi_plus2 + xi_minus2;
      delta = d1/d2;
      delta_mdfy = fmax(pow(delta,2/S),0.1);

      // normalization
      x_plus = (x2d[iptr1]-x2d[iptr2])/xi_plus2;
      z_plus = (z2d[iptr1]-z2d[iptr2])/xi_plus2;
      x_minus = (x2d[iptr3]-x2d[iptr2])/xi_minus2;
      z_minus = (z2d[iptr3]-z2d[iptr2])/xi_minus2;

      dot = x_plus*x_minus + z_plus*z_minus;
      det = x_plus*z_minus - z_plus*x_minus;

      // cal two normal vector clockwise angle. 
      // the method from website
      // from right vector to left vector
      // z axis upward, so is -det
      theta = atan2(-det,dot);
      if(theta<0)
      {
        theta = theta + 2*PI;
      }
      if(theta<=PI)
      {
        alpha = 1.0/(1-pow(cos(theta/2),2));
      }
      if(theta>PI)
      {
        alpha = 1;
      }
      coef_e[i-1] = coef*N_xi*S*delta_mdfy*alpha;
    }
  }

  return 0;
}

int 
cal_matrix(float *x2d, float *z2d, int nx, int k, float *step,
           double *a, double *b, double *c, double *d, float *area)
{
  double A[2][2], B[2][2];
  double mat[2][2], vec[2];
  double mat_b[2][2], vec_d[2];
  size_t iptr1,iptr2,iptr3;
  float x_xi0,z_xi0,x_zt0,z_zt0;
  float x_plus,z_plus,x_minus,z_minus;
  float arc_plus,arc_minus,arc_len,temp;

  if(k==1)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i+1;
      iptr2 = (k-1)*nx + i;
      iptr3 = (k-1)*nx + i-1;
      x_xi0 = 0.5*(x2d[iptr1] - x2d[iptr3]);
      z_xi0 = 0.5*(z2d[iptr1] - z2d[iptr3]);
      x_plus = x2d[iptr1] - x2d[iptr2];
      z_plus = z2d[iptr1] - z2d[iptr2];
      arc_plus = sqrt(pow(x_plus,2) + pow(z_plus,2));
      x_minus = x2d[iptr3] - x2d[iptr2];
      z_minus = z2d[iptr3] - z2d[iptr2];
      arc_minus = sqrt(pow(x_minus,2) + pow(z_minus,2));
      arc_len = 0.5*(arc_plus + arc_minus);
      // arc_length -> area
      // area(i) = A1 
      // area(i+nx) = A0
      area[i] = arc_len * step[k-1];
      area[i+nx] = area[i];
      temp = pow(x_xi0,2) + pow(z_xi0,2);
      x_zt0 = -z_xi0*area[i+nx]/temp;
      z_zt0 =  x_xi0*area[i+nx]/temp;
      A[0][0] = x_zt0; A[0][1] = z_zt0;
      A[1][0] = z_zt0; A[1][1] =-x_zt0; 
      B[0][0] = x_xi0; B[0][1] = z_xi0;
      B[1][0] =-z_xi0; B[1][1] = x_xi0; 
      mat_invert2x2(B);
      mat_mul2x2(B,A,mat);
      mat_iden2x2(mat_b);
      vec[0] = 0; vec[1] = area[i];
      mat_mul2x1(B,vec,vec_d);
      iptr1 = (i-1)*CONST_NDIM*CONST_NDIM;
      iptr2 = (i-1)*CONST_NDIM;
      for(int ii=0; ii<2; ii++) {
        for(int jj=0; jj<2; jj++) {
          a[iptr1+2*ii+jj] = -0.5*mat[ii][jj];
          b[iptr1+2*ii+jj] = mat_b[ii][jj];
          c[iptr1+2*ii+jj] = 0.5*mat[ii][jj];
        }
        d[iptr2+ii] = vec_d[ii];
      }
    }
  } else {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i+1;
      iptr2 = (k-1)*nx + i;
      iptr3 = (k-1)*nx + i-1;
      x_xi0 = 0.5*(x2d[iptr1] - x2d[iptr3]);
      z_xi0 = 0.5*(z2d[iptr1] - z2d[iptr3]);
      x_plus = x2d[iptr1] - x2d[iptr2];
      z_plus = z2d[iptr1] - z2d[iptr2];
      arc_plus = sqrt(pow(x_plus,2) + pow(z_plus,2));
      x_minus = x2d[iptr3] - x2d[iptr2];
      z_minus = z2d[iptr3] - z2d[iptr2];
      arc_minus = sqrt(pow(x_minus,2) + pow(z_minus,2));
      arc_len = 0.5*(arc_plus + arc_minus);
      // arc_length -> area
      // area(i) = A1 
      // area(i+nx) = A0
      area[i+nx] = area[i];
      area[i] = arc_len * step[k-1];
      temp = pow(x_xi0,2) + pow(z_xi0,2);
      x_zt0 = -z_xi0*area[i+nx]/temp;
      z_zt0 =  x_xi0*area[i+nx]/temp;
      A[0][0] = x_zt0; A[0][1] = z_zt0;
      A[1][0] = z_zt0; A[1][1] =-x_zt0; 
      B[0][0] = x_xi0; B[0][1] = z_xi0;
      B[1][0] =-z_xi0; B[1][1] = x_xi0; 
      mat_invert2x2(B);
      mat_mul2x2(B,A,mat);
      mat_iden2x2(mat_b);
      vec[0] = 0; vec[1] = area[i];
      mat_mul2x1(B,vec,vec_d);
      iptr1 = (i-1)*CONST_NDIM*CONST_NDIM;
      iptr2 = (i-1)*CONST_NDIM;
      for(int ii=0; ii<2; ii++) {
        for(int jj=0; jj<2; jj++) {
          a[iptr1+2*ii+jj] = -0.5*mat[ii][jj];
          b[iptr1+2*ii+jj] = mat_b[ii][jj];
          c[iptr1+2*ii+jj] = 0.5*mat[ii][jj];
        }
        d[iptr2+ii] = vec_d[ii];
      }
    }
  }

  return 0;
}

int
modify_smooth(float *x2d, float *z2d, int nx, int k, double *a,
              double *b, double *c, double *d, float *coef_e)
{
  double mat[2][2], vec1[2], vec2[2], vec3[2];
  float coef_i;
  size_t iptr1, iptr2, iptr3, iptr4, iptr5;
  mat_iden2x2(mat);
  for(int i=1; i<nx-1; i++)
  {
    iptr1 = (k-1)*nx + i-1;
    iptr2 = (k-1)*nx + i;
    iptr3 = (k-1)*nx + i+1;
    vec1[0] = x2d[iptr1];
    vec1[1] = z2d[iptr1];
    vec2[0] = x2d[iptr2];
    vec2[1] = z2d[iptr2];
    vec3[0] = x2d[iptr3];
    vec3[1] = z2d[iptr3];

    iptr4 = (i-1)*CONST_NDIM*CONST_NDIM;
    iptr5 = (i-1)*CONST_NDIM;
    
    coef_i = 2*coef_e[i-1];

    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        a[iptr4+2*ii+jj] = a[iptr4+2*ii+jj] - coef_i*mat[ii][jj];
        b[iptr4+2*ii+jj] = b[iptr4+2*ii+jj] + 2*coef_i*mat[ii][jj];
        c[iptr4+2*ii+jj] = c[iptr4+2*ii+jj] - coef_i*mat[ii][jj];
      }
      d[iptr5+ii] = d[iptr5+ii] + coef_e[i-1]*(vec1[ii]+vec3[ii]-2*vec2[ii]);
    }
  }

  return 0;
}

int
modify_bdry(int n, double *a, double *b, double *c, double *d,
            float epsilon, int bdry_itype)
{
  size_t iptr;
  // float boundry
  if(bdry_itype == 1)
  {
    iptr = (n-1)*CONST_NDIM*CONST_NDIM;
    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        // modify i=0
        b[ii*2+jj] = b[ii*2+jj] + (1+epsilon)*a[ii*2+jj];
        c[ii*2+jj] = c[ii*2+jj] - epsilon*a[ii*2+jj];
        // modify i=n-1
        b[iptr+ii*2+jj] = b[iptr+ii*2+jj] + (1+epsilon)*c[iptr+ii*2+jj];
        a[iptr+ii*2+jj] = a[iptr+ii*2+jj] - epsilon*c[iptr+ii*2+jj];
      }
    }
  }

  // dx=0 cartesian boundry
  if(bdry_itype == 2)
  {
    iptr = (n-1)*CONST_NDIM*CONST_NDIM;
    for(int ii=0; ii<2; ii++) {
      // only modify second column
      int jj = 1;
      // modify i=0
      b[ii*2+jj] = b[ii*2+jj] + a[ii*2+jj];
      // modify i=nx-3
      b[iptr+ii*2+jj] = b[iptr+ii*2+jj] + c[iptr+ii*2+jj];
    }
  }

  return 0;
}

int
assign_coords(double *xz, float *x2d, float *z2d, int nx, int k,
              float epsilon, int bdry_itype)
{
  size_t iptr,iptr1,iptr2;
  size_t iptr3,iptr4,iptr5;
  for(int i=1; i<nx-1; i++)
  {
    iptr  =  k*nx + i;
    iptr1 = (k-1)*nx + i;
    iptr2 = (i-1)*CONST_NDIM;
    x2d[iptr] = x2d[iptr1] + xz[iptr2];
    z2d[iptr] = z2d[iptr1] + xz[iptr2+1];
  }

  // floating boundary
  if(bdry_itype == 1)
  {
    iptr  = k*nx+0;       // (0,k)
    iptr1 = (k-1)*nx+0;   // (0,k-1)
    iptr2 = k*nx+1;       // (1,k)
    iptr3 = (k-1)*nx+1;   // (1,k-1)
    iptr4 = k*nx+2;       // (2,k)
    iptr5 = (k-1)*nx+2;   // (2,k-1)

    x2d[iptr] = x2d[iptr1] + (1+epsilon)*(x2d[iptr2]-x2d[iptr3])
               -epsilon*(x2d[iptr4]-x2d[iptr5]);
    z2d[iptr] = z2d[iptr1] + (1+epsilon)*(z2d[iptr2]-z2d[iptr3])
               -epsilon*(z2d[iptr4]-z2d[iptr5]);

    iptr  = k*nx+(nx-1);       // (nx-1,k)
    iptr1 = (k-1)*nx+(nx-1);   // (nx-1,k-1)
    iptr2 = k*nx+(nx-2);       // (nx-2,k)
    iptr3 = (k-1)*nx+(nx-2);   // (nx-2,k-1)
    iptr4 = k*nx+(nx-3);       // (nx-3,k)
    iptr5 = (k-1)*nx+(nx-3);   // (nx-3,k-1)

    x2d[iptr] = x2d[iptr1] + (1+epsilon)*(x2d[iptr2]-x2d[iptr3])
               -epsilon*(x2d[iptr4]-x2d[iptr5]);
    z2d[iptr] = z2d[iptr1] + (1+epsilon)*(z2d[iptr2]-z2d[iptr3])
               -epsilon*(z2d[iptr4]-z2d[iptr5]);
  }

  // cartesian boundary
  if(bdry_itype == 2)
  {
    iptr = k*nx+0;        // (0,k) 
    iptr1 = (k-1)*nx+0;   // (0,k-1)
    iptr2 = k*nx+1;       // (1,k)
    iptr3 = (k-1)*nx+1;   // (1,k-1)

    x2d[iptr] = x2d[iptr1];
    z2d[iptr] = z2d[iptr1] + z2d[iptr2] - z2d[iptr3];

    iptr = k*nx+(nx-1);        // (nx-1,k) 
    iptr1 = (k-1)*nx+(nx-1);   // (nx-1,k-1)
    iptr2 = k*nx+1;            // (nx-2,k)
    iptr3 = (k-1)*nx+1;        // (nx-2,k-1)

    x2d[iptr] = x2d[iptr1];
    z2d[iptr] = z2d[iptr1] + z2d[iptr2] - z2d[iptr3];
  }

  return 0;
}

