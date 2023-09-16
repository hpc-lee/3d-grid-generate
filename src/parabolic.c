#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "parabolic.h"
#include "solver.h"
#include "lib_mem.h"

int 
para_gene(gd_t *gdcurv, float coef, int o2i)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  if(o2i == 1)
  {
    fprintf(stdout,"parabolic method, march direction from inner bdry(k=0) to outer bdry(nz-1)\n");
    fprintf(stdout,"so if we want from outer bdry(nz-1) to inner bdry(k=0), must be flip\n");
    fprintf(stdout,"flip coord, trans outer bdry(nz-1) to inner bdry(0)\n");
    flip_coord_z(gdcurv);
  }
  // malloc space for thomas method 
  // a,b,c,d_x,d_y,d_z,u_x,u_y,u_z. 9 vars space
  float *var_th = NULL;
  var_th = (float *)mem_calloc_1d_float((nx-2)*9, 0.0, "init");

  // calculate x1 x2 y1 y2 bdry arc length
  float *x1_len = NULL;
  x1_len = (float *)mem_calloc_1d_float(nz*ny, 0.0, "init");
  float *x2_len = NULL;
  x2_len = (float *)mem_calloc_1d_float(nz*ny, 0.0, "init");
  float *y1_len = NULL;
  y1_len = (float *)mem_calloc_1d_float(nz*nx, 0.0, "init");
  float *y2_len = NULL;
  y2_len = (float *)mem_calloc_1d_float(nz*nx, 0.0, "init");

  cal_bdry_arc_length(x3d,y3d,z3d,nx,ny,nz,x1_len,x2_len,y1_len,y2_len);

  for(int k=1; k<nz-1; k++)
  {
    // predict k+1 layer points
    predict_point(x3d,y3d,z3d,nx,ny,nz,k,o2i,coef,x1_len,x2_len,y1_len,y2_len);
    // base predict points
    // update k layer points
    update_point(x3d,y3d,z3d,var_th,nx,ny,nz,k);
    fprintf(stdout,"number of layer is %d\n",k);
    fflush(stdout);
  }

  if(o2i == 1)
  {
    // flip to return.
    flip_coord_z(gdcurv);
  }

  free(var_th);
  free(x1_len);
  free(x2_len);
  free(y1_len);
  free(y2_len);

  return 0;
}

int 
predict_point(float *x3d, float *y3d, float *z3d, int nx, int ny, 
              int nz, int k, int o2i, float coef, float *x1_len,
              float *x2_len, float *y1_len, float *y2_len)
{
  // k-1 layer points is know
  // predict points k+1 and k layer
  // first predict further k+1 layer
  // base predict k+1 layer, calculate k layer
  // NOTE: this predict point only used
  // by calculate matrix or coefficient
  // not the final point
  //
  // cal k-1 layer point unit normal vector
  // vn -> vector normal
  
  float xi,et,zt,cs;
  size_t iptr1,iptr2,iptr3,iptr4;
  float vn_x,vn_y,vn_z,vn_len;
  float vn_x0,vn_y0,vn_z0;
  float R_x,R_y,R_z,R;
  float R_x1, r_x1, r_x11;
  float R_x2, r_x2, r_x22;
  float R_y1, r_y1, r_y11;
  float R_y2, r_y2, r_y22;
  float c,cc,c_x,cc_x,c_y,cc_y;
  float c_x1,c_x11,c_x2,c_x22;
  float c_y1,c_y11,c_y2,c_y22;
  float x0,y0,z0,xs,ys,zs;
  float x_xi, y_xi, z_xi;
  float x_et, y_et, z_et;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;

  zt = (1.0*(k-1))/((nz-1)-2);
  cs = exp(coef*zt);

  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      // cal normal vector 
      iptr1 = (k-1)*siz_iz + j*siz_iy + i+1;
      iptr2 = (k-1)*siz_iz + j*siz_iy + i-1;
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);

      iptr3 = (k-1)*siz_iz + (j+1)*siz_iy + i;
      iptr4 = (k-1)*siz_iz + (j-1)*siz_iy + i;
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);

      vn_x = y_xi*z_et - z_xi*y_et; 
      vn_y = z_xi*x_et - x_xi*z_et; 
      vn_z = x_xi*y_et - y_xi*x_et; 
      if(o2i == 1)
      {
        vn_x = -vn_x;
        vn_y = -vn_y;
        vn_z = -vn_z;
      } 
      vn_len = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      vn_x0 = vn_x/vn_len;
      vn_y0 = vn_y/vn_len;
      vn_z0 = vn_z/vn_len;
      
      // inner point i
      iptr1 = (k-1)*siz_iz + j*siz_iy + i; // (i,j,k-1)
      iptr2 = (nz-1)*siz_iz + j*siz_iy + i; // (i,j,nz-1)
      R_x = x3d[iptr1] - x3d[iptr2];
      R_y = y3d[iptr1] - y3d[iptr2];
      R_z = z3d[iptr1] - z3d[iptr2];
      R = sqrt(pow(R_x,2)+pow(R_y,2)+pow(R_z,2));

      // bdry point i=0
      iptr1 = (k-1)*ny + j;
      iptr2 = (nz-1)*ny + j;
      iptr3 = (k+1)*ny + j;
      iptr4 = k*ny + j;
      R_x1  = x1_len[iptr2] - x1_len[iptr1];
      r_x1  = x1_len[iptr3] - x1_len[iptr1];
      r_x11 = x1_len[iptr4] - x1_len[iptr1];

      c_x1  = r_x1/R_x1;
      c_x11 = r_x11/r_x1;

      // bdry point i=nx-1
      iptr1 = (k-1)*ny + j;
      iptr2 = (nz-1)*ny + j;
      iptr3 = (k+1)*ny + j;
      iptr4 = k*ny + j;
      R_x2  = x2_len[iptr2] - x2_len[iptr1];
      r_x2  = x2_len[iptr3] - x2_len[iptr1];
      r_x22 = x2_len[iptr4] - x2_len[iptr1];

      c_x2  = r_x2/R_x2;
      c_x22 = r_x22/r_x2;

      // cal clustering factor c
      xi = (1.0*i)/(nx-1);
      c_x =  (1-xi)*c_x1  + xi*c_x2;
      cc_x = (1-xi)*c_x11 + xi*c_x22;
      
      // bdry point j=0
      iptr1 = (k-1)*nx + i;
      iptr2 = (nz-1)*nx + i;
      iptr3 = (k+1)*nx + i;
      iptr4 = k*nx + i;
      R_y1  = y1_len[iptr2] - y1_len[iptr1];
      r_y1  = y1_len[iptr3] - y1_len[iptr1];
      r_y11 = y1_len[iptr4] - y1_len[iptr1];

      c_y1  = r_y1/R_y1;
      c_y11 = r_y11/r_y1;

      // bdry point j=ny-1
      iptr1 = (k-1)*nx + i;
      iptr2 = (nz-1)*nx + i;
      iptr3 = (k+1)*nx + i;
      iptr4 = k*nx + i;
      R_y2  = y2_len[iptr2] - y2_len[iptr1];
      r_y2  = y2_len[iptr3] - y2_len[iptr1];
      r_y22 = y2_len[iptr4] - y2_len[iptr1];

      c_y2  = r_y2/R_y2;
      c_y22 = r_y22/r_y2;

      // cal clustering factor c
      et = (1.0*j)/(ny-1);
      c_y =  (1-et)*c_y1  + et*c_y2;
      cc_y = (1-et)*c_y11 + et*c_y22;

      c = 0.5*(c_x+c_y);
      cc = 0.5*(cc_x+cc_y);

      iptr1 = (k-1)*siz_iz + j*siz_iy + i; //(i,j,k-1)
      // x0,y0,z0 normal point
      x0 = x3d[iptr1] + vn_x0*c*R;
      y0 = y3d[iptr1] + vn_y0*c*R;
      z0 = z3d[iptr1] + vn_z0*c*R;
      
      // xs,ys,zs linear distance point
      iptr2 = (nz-1)*siz_iz + j*siz_iy + i;   //(i,j,nz-1)
      xs = x3d[iptr1] + c*(x3d[iptr2]-x3d[iptr1]);
      ys = y3d[iptr1] + c*(y3d[iptr2]-y3d[iptr1]);
      zs = z3d[iptr1] + c*(z3d[iptr2]-z3d[iptr1]);

      iptr3 = (k+1)*siz_iz + j*siz_iy + i;     //(i,j,k+1)
      x3d[iptr3] = cs*x0 + (1-cs)*xs;
      y3d[iptr3] = cs*y0 + (1-cs)*ys;
      z3d[iptr3] = cs*z0 + (1-cs)*zs;

      iptr4 = k*siz_iz + j*siz_iy + i;    //(i,j,k)
      x3d[iptr4] = x3d[iptr1] + cc*(x3d[iptr3]-x3d[iptr1]);
      y3d[iptr4] = y3d[iptr1] + cc*(y3d[iptr3]-y3d[iptr1]);
      z3d[iptr4] = z3d[iptr1] + cc*(z3d[iptr3]-z3d[iptr1]);
    }
  }

  //bdry_effct(x3d,y3d,z3d,nx,ny,nz,k);

  return 0;
}

int
update_point(float *x3d, float *y3d, float *z3d, float *var_th, int nx, int ny, int nz, int k)
{
  float *a;
  float *b;
  float *c;
  float *d_x;
  float *d_y;
  float *d_z;
  float *u_x;
  float *u_y;
  float *u_z;
  int siz_vec = nx-2;
  a = var_th;
  b = var_th + siz_vec;
  c = var_th + 2*siz_vec;
  d_x = var_th + 3*siz_vec;
  d_y = var_th + 4*siz_vec;
  d_z = var_th + 5*siz_vec;
  u_x = var_th + 6*siz_vec;
  u_y = var_th + 7*siz_vec;
  u_z = var_th + 8*siz_vec;

  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  size_t iptr5,iptr6;
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xiet,y_xiet,z_xiet;
  float x_xizt,y_xizt,z_xizt;
  float x_etzt,y_etzt,z_etzt;
  float g11,g22,g33,g12,g13,g23;
  float alpha1,alpha2,alpha3;
  float beta12,beta23,beta13;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;

  for(int j=1; j<ny-1; j++) 
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = k*siz_iz + j*siz_iy + (i+1);  // (i+1,j,k)
      iptr2 = k*siz_iz + j*siz_iy + (i-1);  // (i-1,j,k)
      iptr3 = k*siz_iz + (j+1)*siz_iy + i;  // (i,j+1,k)
      iptr4 = k*siz_iz + (j-1)*siz_iy + i;  // (i,j-1,k)
      iptr5 = (k+1)*siz_iz + j*siz_iy + i;  // (i,j,k+1)
      iptr6 = (k-1)*siz_iz + j*siz_iy + i;  // (i,j,k-1)

      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);
      x_zt = 0.5*(x3d[iptr5] - x3d[iptr6]);
      y_zt = 0.5*(y3d[iptr5] - y3d[iptr6]);
      z_zt = 0.5*(z3d[iptr5] - z3d[iptr6]);

      iptr1 = k*siz_iz + (j-1)*siz_iy + (i-1);   // (i-1,j-1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + (i+1);   // (i+1,j-1,k)
      iptr3 = k*siz_iz + (j+1)*siz_iy + (i+1);   // (i+1,j+1,k)
      iptr4 = k*siz_iz + (j+1)*siz_iy + (i-1);   // (i-1,j+1,k)

      x_xiet = 0.25*(x3d[iptr3] + x3d[iptr1] - x3d[iptr2] - x3d[iptr4]);
      y_xiet = 0.25*(y3d[iptr3] + y3d[iptr1] - y3d[iptr2] - y3d[iptr4]);
      z_xiet = 0.25*(z3d[iptr3] + z3d[iptr1] - z3d[iptr2] - z3d[iptr4]);

      iptr1 = (k-1)*siz_iz + j*siz_iy + (i-1);   // (i-1,j,k-1)
      iptr2 = (k-1)*siz_iz + j*siz_iy + (i+1);   // (i+1,j,k-1)
      iptr3 = (k+1)*siz_iz + j*siz_iy + (i+1);   // (i+1,j,k+1)
      iptr4 = (k+1)*siz_iz + j*siz_iy + (i-1);   // (i-1,j,k+1)

      x_xizt = 0.25*(x3d[iptr3] + x3d[iptr1] - x3d[iptr2] - x3d[iptr4]);
      y_xizt = 0.25*(y3d[iptr3] + y3d[iptr1] - y3d[iptr2] - y3d[iptr4]);
      z_xizt = 0.25*(z3d[iptr3] + z3d[iptr1] - z3d[iptr2] - z3d[iptr4]);

      iptr1 = (k-1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-1)
      iptr2 = (k-1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-1)
      iptr3 = (k+1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k+1)
      iptr4 = (k+1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k+1)

      x_etzt = 0.25*(x3d[iptr3] + x3d[iptr1] - x3d[iptr2] - x3d[iptr4]);
      y_etzt = 0.25*(y3d[iptr3] + y3d[iptr1] - y3d[iptr2] - y3d[iptr4]);
      z_etzt = 0.25*(z3d[iptr3] + z3d[iptr1] - z3d[iptr2] - z3d[iptr4]);

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

      a[i-1] = alpha1;
      b[i-1] = -2*(alpha1 + alpha2 + alpha3);
      c[i-1] = alpha1;

      iptr1 = k*siz_iz + (j+1)*siz_iy + i;
      iptr2 = k*siz_iz + (j-1)*siz_iy + i;
      iptr3 = (k+1)*siz_iz + j*siz_iy + i;
      iptr4 = (k-1)*siz_iz + j*siz_iy + i;

      d_x[i-1] = -alpha2*(x3d[iptr1]+x3d[iptr2]) - alpha3*(x3d[iptr3]+x3d[iptr4]) 
                 -2*(beta12*x_xiet + beta23*x_etzt + beta13*x_xizt);

      d_y[i-1] = -alpha2*(y3d[iptr1]+y3d[iptr2]) - alpha3*(y3d[iptr3]+y3d[iptr4]) 
                 -2*(beta12*y_xiet + beta23*y_etzt + beta13*y_xizt);

      d_z[i-1] = -alpha2*(z3d[iptr1]+z3d[iptr2]) - alpha3*(z3d[iptr3]+z3d[iptr4]) 
                 -2*(beta12*z_xiet + beta23*z_etzt + beta13*z_xizt);
    }

    // i=1 modify a,d_x,d_y,d_z
    iptr = k*siz_iz + j*siz_iy + 0;  // (0,j,k)
    d_x[0] = d_x[0] - a[0]*x3d[iptr];
    d_y[0] = d_y[0] - a[0]*y3d[iptr];
    d_z[0] = d_z[0] - a[0]*z3d[iptr];
    a[0] = 0;

    // i=nx-2 modify c,d_x,d_z
    iptr = k*siz_iz + j*siz_iy + (nx-1);
    d_x[nx-3] = d_x[nx-3] - c[nx-3]*x3d[iptr];
    d_y[nx-3] = d_y[nx-3] - c[nx-3]*y3d[iptr];
    d_z[nx-3] = d_z[nx-3] - c[nx-3]*z3d[iptr];
    c[nx-3] = 0;
    
    // cal coords and update
    thomas(siz_vec,a,b,c,d_x,d_y,d_z,u_x,u_y,u_z);
    for(int i=1; i<nx-1; i++)
    {
      iptr = k*siz_iz + j*siz_iy + i;
      x3d[iptr] = u_x[i-1];
      y3d[iptr] = u_y[i-1];
      z3d[iptr] = u_z[i-1];
    }
  }

  return 0;
}
/*
int
bdry_effct(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, int k)
{
  // add bdry point effct
  // due to bdry maybe nonlinear,
  // bdry undulate
  float b1_x1_k, b1_z1_k, b1_x1_k1, b1_z1_k1;
  float b1_x2_k, b1_z2_k, b1_x2_k1, b1_z2_k1;
  float b2_x1_k, b2_z1_k, b2_x1_k1, b2_z1_k1;
  float b2_x2_k, b2_z2_k, b2_x2_k1, b2_z2_k1;
  float dif_b1_x_k, dif_b1_z_k, dif_b1_x_k1, dif_b1_z_k1;
  float dif_b2_x_k, dif_b2_z_k, dif_b2_x_k1, dif_b2_z_k1;
  size_t iptr1,iptr2,iptr3;
  // x1 bdry
  iptr1 = (k-1)*nx + 0; //(0,k-1)
  iptr2 =  k*nx + 0;    //(0,k)
  iptr3 = (k+1)*nx + 0; //(0,k+1)
  b1_x1_k  = x3d[iptr2]-x3d[iptr1];
  b1_z1_k  = z3d[iptr2]-z3d[iptr1];
  b1_x1_k1 = x3d[iptr3]-x3d[iptr2];
  b1_z1_k1 = z3d[iptr3]-z3d[iptr2];

  iptr1 = (k-1)*nx + 1; //(1,k-1)
  iptr2 =  k*nx + 1;    //(1,k)
  iptr3 = (k+1)*nx + 1; //(1,k+1)
  b1_x2_k  = x3d[iptr2]-x3d[iptr1];
  b1_z2_k  = z3d[iptr2]-z3d[iptr1];
  b1_x2_k1 = x3d[iptr3]-x3d[iptr2];
  b1_z2_k1 = z3d[iptr3]-z3d[iptr2];

  dif_b1_x_k  = b1_x1_k -b1_x2_k;
  dif_b1_z_k  = b1_z1_k -b1_z2_k;
  dif_b1_x_k1 = b1_x1_k1-b1_x2_k1;
  dif_b1_z_k1 = b1_z1_k1-b1_z2_k1;

  // x2 bdry
  iptr1 = (k-1)*nx + (nx-1); //(nx-1,k-1)
  iptr2 =  k*nx + (nx-1);    //(nx-1,k)
  iptr3 = (k+1)*nx + (nx-1); //(nx-1,k+1)
  b2_x1_k  = x3d[iptr2]-x3d[iptr1];
  b2_z1_k  = z3d[iptr2]-z3d[iptr1];
  b2_x1_k1 = x3d[iptr3]-x3d[iptr2];
  b2_z1_k1 = z3d[iptr3]-z3d[iptr2];

  iptr1 = (k-1)*nx + (nx-2); //(nx-2,k-1)
  iptr2 =  k*nx + (nx-2);    //(nx-2,k)
  iptr3 = (k+1)*nx + (nx-2); //(nx-2,k+1)
  b2_x2_k  = x3d[iptr2]-x3d[iptr1];
  b2_z2_k  = z3d[iptr2]-z3d[iptr1];
  b2_x2_k1 = x3d[iptr3]-x3d[iptr2];
  b2_z2_k1 = z3d[iptr3]-z3d[iptr2];
  
  dif_b2_x_k  = b2_x1_k -b2_x2_k;
  dif_b2_z_k  = b2_z1_k -b2_z2_k;
  dif_b2_x_k1 = b2_x1_k1-b2_x2_k1;
  dif_b2_z_k1 = b2_z1_k1-b2_z2_k1;

  float xi, bdry_x_k, bdry_z_k;
  float bdry_x_k1, bdry_z_k1;
  for(int i=1; i<nx-1; i++)
  {
    xi = (1.0*i)/(nx-1);
    bdry_x_k = (1-xi)*dif_b1_x_k + xi*dif_b2_x_k;
    bdry_z_k = (1-xi)*dif_b1_z_k + xi*dif_b2_z_k;
    iptr1 = k*nx + i; 
    x3d[iptr1] = x3d[iptr1] + 2.05*bdry_x_k;
    z3d[iptr1] = z3d[iptr1] + 2.05*bdry_z_k;

    bdry_x_k1 = (1-xi)*dif_b1_x_k1 + xi*dif_b2_x_k1;
    bdry_z_k1 = (1-xi)*dif_b1_z_k1 + xi*dif_b2_z_k1;
    iptr2 = (k+1)*nx + i; 
    x3d[iptr2] = x3d[iptr2] + 2.05*bdry_x_k1;
    z3d[iptr2] = z3d[iptr2] + 2.05*bdry_z_k1;
  }

  return 0;
}
*/

int
cal_bdry_arc_length(float *x3d, float *y3d, float *z3d, int nx,
                    int ny, int nz, float *x1_len, float *x2_len,
                    float *y1_len, float *y2_len)
{
  size_t iptr1, iptr2, iptr3, iptr4;
  float x_len, y_len, z_len, dh_len;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;
  for(int k=1; k<nz; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr1 = k*siz_iz + j*siz_iy + 0;
      iptr2 = (k-1)*siz_iz + j*siz_iy + 0;
      x_len = x3d[iptr1] - x3d[iptr2];
      y_len = y3d[iptr1] - y3d[iptr2];
      z_len = z3d[iptr1] - z3d[iptr2];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      iptr3 = k*ny + j;
      iptr4 = (k-1)*ny + j;
      x1_len[iptr3] = x1_len[iptr4] + dh_len;

      iptr1 = k*siz_iz + j*siz_iy + nx-1;
      iptr2 = (k-1)*siz_iz + j*siz_iy + nx-1;
      x_len = x3d[iptr1] - x3d[iptr2];
      y_len = y3d[iptr1] - y3d[iptr2];
      z_len = z3d[iptr1] - z3d[iptr2];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      iptr3 = k*ny + j;
      iptr4 = (k-1)*ny + j;
      x2_len[iptr3] = x2_len[iptr4] + dh_len;
    }
  }

  for(int k=1; k<nz; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = k*siz_iz + 0*siz_iy + i;
      iptr2 = (k-1)*siz_iz + 0*siz_iy + i;
      x_len = x3d[iptr1] - x3d[iptr2];
      y_len = y3d[iptr1] - y3d[iptr2];
      z_len = z3d[iptr1] - z3d[iptr2];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      iptr3 = k*nx + i;
      iptr4 = (k-1)*nx + i;
      y1_len[iptr3] = y1_len[iptr4] + dh_len;

      iptr1 = k*siz_iz + (ny-1)*siz_iy + i;
      iptr2 = (k-1)*siz_iz + (ny-1)*siz_iy + i;
      x_len = x3d[iptr1] - x3d[iptr2];
      y_len = y3d[iptr1] - y3d[iptr2];
      z_len = z3d[iptr1] - z3d[iptr2];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      iptr3 = k*nx + i;
      iptr4 = (k-1)*nx + i;
      y2_len[iptr3] = y2_len[iptr4] + dh_len;
    }
  }
  return 0;
}
