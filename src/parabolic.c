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
  int nz = gdcurv->nz;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  if(o2i == 1)
  {
    fprintf(stdout,"parabolic method, march direction from inner bdry(k=0) to outer bdry(nz-1)\n");
    fprintf(stdout,"so if we want from outer bdry(nz-1) to inner bdry(k=0), must be flip\n");
    fprintf(stdout,"flip coord, trans outer bdry(nz-1) to inner bdry(0)\n");
    flip_coord(x2d,nx,nz);
    flip_coord(z2d,nx,nz);
  }
  // malloc space for thomas method 
  // a,b,c,d_x,d_z,u_x,u_z. 7 vars space
  float *var_th = NULL;
  var_th = (float *)mem_calloc_1d_float((nx-2)*7, 0.0, "init");

  // calculate x1 x2 bdry arc length
  float *x1_len = NULL;
  x1_len = (float *)mem_calloc_1d_float(nz, 0.0, "init");
  float *x2_len = NULL;
  x2_len = (float *)mem_calloc_1d_float(nz, 0.0, "init");

  cal_bdry_arc_length(x2d,z2d,nx,nz,x1_len,x2_len);

  for(int k=1; k<nz-1; k++)
  {
    // predict k+1 layer points
    predict_point(x2d,z2d,nx,nz,k,o2i,coef,x1_len,x2_len);
    // base predict points
    // update k layer points
    update_point(x2d,z2d,var_th,nx,nz,k);
  }

  if(o2i == 1)
  {
    // flip to return.
    flip_coord(x2d,nx,nz);
    flip_coord(z2d,nx,nz);
  }

  free(var_th);
  free(x1_len);
  free(x2_len);
  return 0;
}

int 
predict_point(float *x2d, float *z2d, int nx, int nz, int k, int o2i, 
              float coef, float *x1_len, float *x2_len)
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
  // vt -> vector tangential
  // vn -> vector normal
  
  float xi,zt,cs;
  size_t iptr1,iptr2,iptr3,iptr4;
  int sign1;
  float vt_x,vt_z,vn_x,vn_z,len_vt;
  float R_x,R_z,R;
  float R1, r1, r11;
  float R2, r2, r22;
  float c,cc,c1,c11,c2,c22;
  float x0,z0,xs,zs;

  zt = (1.0*(k-1))/((nz-1)-2);
  cs = exp(coef*zt);
  // o2i outer bdry to inner bdry
  if(o2i == 1) { 
    sign1 = 1;
  } else {
    sign1 = -1;
  }

  for(int i=1; i<nx-1; i++)
  {
    // cal normal vector 
    iptr1 = (k-1)*nx + (i+1); // (i+1,k-1)
    iptr2 = (k-1)*nx + (i-1); // (i-1,k-1)
    vt_x = 0.5*(x2d[iptr1] - x2d[iptr2]); 
    vt_z = 0.5*(z2d[iptr1] - z2d[iptr2]); 
    len_vt = sqrt(pow(vt_x,2)+pow(vt_z,2));
    vn_x =  sign1*vt_z/len_vt;
    vn_z = -sign1*vt_x/len_vt;
    
    // inner point i
    iptr1 = (k-1)*nx + i; // (i,k-1)
    iptr2 = (nz-1)*nx + i; // (i,nz-1)
    R_x = x2d[iptr1] - x2d[iptr2];
    R_z = z2d[iptr1] - z2d[iptr2];
    R = sqrt(pow(R_x,2)+pow(R_z,2));

    // bdry point i=0
    R1 = x1_len[nz-1] - x1_len[k-1];

    r1 = x1_len[k+1] - x1_len[k-1];

    r11 = x1_len[k] - x1_len[k-1];

    c1 = r1/R1;
    c11 = r11/r1;

    // bdry point i=nx-1
    R2 = x2_len[nz-1] - x2_len[k-1];

    r2 = x2_len[k+1] - x2_len[k-1];

    r22 = x2_len[k] - x2_len[k-1];

    c2 = r2/R2;
    c22 = r22/r2;

    // cal clustering factor c
    xi = (1.0*i)/(nx-1);
    c =  (1-xi)*c1  + xi*c2;
    cc = (1-xi)*c11 + xi*c22;
    
    iptr1 = (k-1)*nx + i; //(i,k-1)
    // x0,z0 normal point
    x0 = x2d[iptr1] + vn_x*c*R;
    z0 = z2d[iptr1] + vn_z*c*R;
    
    // xs,zs linear distance point
    iptr1 = (k-1)*nx + i; //(i,k-1)
    iptr2 = (nz-1)*nx + i;   //(i,nz-1)
    xs = x2d[iptr1] + c*(x2d[iptr2]-x2d[iptr1]);
    zs = z2d[iptr1] + c*(z2d[iptr2]-z2d[iptr1]);

    iptr3 = (k+1)*nx + i;     //(i,k+1)
    x2d[iptr3] = cs*x0 + (1-cs)*xs;
    z2d[iptr3] = cs*z0 + (1-cs)*zs;

    iptr4 = k*nx + i;    //(i,k)
    x2d[iptr4] = x2d[iptr1] + cc*(x2d[iptr3]-x2d[iptr1]);
    z2d[iptr4] = z2d[iptr1] + cc*(z2d[iptr3]-z2d[iptr1]);
  }

  bdry_effct(x2d,z2d,nx,nz,k);

  return 0;
}

int
update_point(float *x2d, float *z2d, float *var_th, int nx, int nz, int k)
{
  float *a;
  float *b;
  float *c;
  float *d_x;
  float *d_z;
  float *u_x;
  float *u_z;
  int siz_vec = nx-2;
  a = var_th;
  b = var_th + siz_vec;
  c = var_th + 2*siz_vec;
  d_x = var_th + 3*siz_vec;
  d_z = var_th + 4*siz_vec;
  u_x = var_th + 5*siz_vec;
  u_z = var_th + 6*siz_vec;

  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float temp_x,temp_z;
  float x_xizt,z_xizt;
  float g11,g22,g12;

  for(int i=1; i<nx-1; i++)
  {
    iptr1 =  k*nx + (i-1); // (i-1,k)
    iptr2 =  k*nx + (i+1); // (i+1,k)
    iptr3 = (k-1)*nx + i;  // (i,k-1)
    iptr4 = (k+1)*nx + i;  // (i,k+1)

    x_xi = 0.5*(x2d[iptr2] - x2d[iptr1]);
    z_xi = 0.5*(z2d[iptr2] - z2d[iptr1]);

    x_zt = 0.5*(x2d[iptr4] - x2d[iptr3]);
    z_zt = 0.5*(z2d[iptr4] - z2d[iptr3]);

    temp_x = x2d[iptr4] + x2d[iptr3];
    temp_z = z2d[iptr4] + z2d[iptr3];

    iptr1 = (k-1)*nx + (i-1);   // (i-1,k-1)
    iptr2 = (k-1)*nx + (i+1);   // (i+1,k-1)
    iptr3 = (k+1)*nx + (i-1);   // (i-1,k+1)
    iptr4 = (k+1)*nx + (i+1);   // (i+1,k+1)

    x_xizt = 0.25*(x2d[iptr4] + x2d[iptr1] - x2d[iptr2] - x2d[iptr3]);
    z_xizt = 0.25*(z2d[iptr4] + z2d[iptr1] - z2d[iptr2] - z2d[iptr3]);

    g11 = x_xi*x_xi + z_xi*z_xi;
    g22 = x_zt*x_zt + z_zt*z_zt;
    g12 = x_xi*x_zt + z_xi*z_zt;

    a[i-1] = g22;
    b[i-1] = -2*(g22+g11);
    c[i-1] = g22;

    d_x[i-1] = -g11*temp_x + 2*g12*x_xizt;
    d_z[i-1] = -g11*temp_z + 2*g12*z_xizt;
  }

  // i=1 modify a,d_x,d_z
  iptr = k*nx+0;
  d_x[0] = d_x[0] - a[0]*x2d[iptr];
  d_z[0] = d_z[0] - a[0]*z2d[iptr];
  a[0] = 0;

  // i=nx-2 modify c,d_x,d_z
  iptr = k*nx + (nx-1);
  d_x[nx-3] = d_x[nx-3] - c[nx-3]*x2d[iptr];
  d_z[nx-3] = d_z[nx-3] - c[nx-3]*z2d[iptr];
  c[nx-3] = 0;
  
  // cal coords and update
  thomas(siz_vec,a,b,c,d_x,d_z,u_x,u_z);
  for(int i=1; i<nx-1; i++)
  {
    iptr = k*nx + i;
    x2d[iptr] = u_x[i-1];
    z2d[iptr] = u_z[i-1];
  }

  return 0;
}

int
bdry_effct(float *x2d, float *z2d, int nx, int nz, int k)
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
  b1_x1_k  = x2d[iptr2]-x2d[iptr1];
  b1_z1_k  = z2d[iptr2]-z2d[iptr1];
  b1_x1_k1 = x2d[iptr3]-x2d[iptr2];
  b1_z1_k1 = z2d[iptr3]-z2d[iptr2];

  iptr1 = (k-1)*nx + 1; //(1,k-1)
  iptr2 =  k*nx + 1;    //(1,k)
  iptr3 = (k+1)*nx + 1; //(1,k+1)
  b1_x2_k  = x2d[iptr2]-x2d[iptr1];
  b1_z2_k  = z2d[iptr2]-z2d[iptr1];
  b1_x2_k1 = x2d[iptr3]-x2d[iptr2];
  b1_z2_k1 = z2d[iptr3]-z2d[iptr2];

  dif_b1_x_k  = b1_x1_k -b1_x2_k;
  dif_b1_z_k  = b1_z1_k -b1_z2_k;
  dif_b1_x_k1 = b1_x1_k1-b1_x2_k1;
  dif_b1_z_k1 = b1_z1_k1-b1_z2_k1;

  // x2 bdry
  iptr1 = (k-1)*nx + (nx-1); //(nx-1,k-1)
  iptr2 =  k*nx + (nx-1);    //(nx-1,k)
  iptr3 = (k+1)*nx + (nx-1); //(nx-1,k+1)
  b2_x1_k  = x2d[iptr2]-x2d[iptr1];
  b2_z1_k  = z2d[iptr2]-z2d[iptr1];
  b2_x1_k1 = x2d[iptr3]-x2d[iptr2];
  b2_z1_k1 = z2d[iptr3]-z2d[iptr2];

  iptr1 = (k-1)*nx + (nx-2); //(nx-2,k-1)
  iptr2 =  k*nx + (nx-2);    //(nx-2,k)
  iptr3 = (k+1)*nx + (nx-2); //(nx-2,k+1)
  b2_x2_k  = x2d[iptr2]-x2d[iptr1];
  b2_z2_k  = z2d[iptr2]-z2d[iptr1];
  b2_x2_k1 = x2d[iptr3]-x2d[iptr2];
  b2_z2_k1 = z2d[iptr3]-z2d[iptr2];
  
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
    x2d[iptr1] = x2d[iptr1] + 2.05*bdry_x_k;
    z2d[iptr1] = z2d[iptr1] + 2.05*bdry_z_k;

    bdry_x_k1 = (1-xi)*dif_b1_x_k1 + xi*dif_b2_x_k1;
    bdry_z_k1 = (1-xi)*dif_b1_z_k1 + xi*dif_b2_z_k1;
    iptr2 = (k+1)*nx + i; 
    x2d[iptr2] = x2d[iptr2] + 2.05*bdry_x_k1;
    z2d[iptr2] = z2d[iptr2] + 2.05*bdry_z_k1;
  }

  return 0;
}

int
cal_bdry_arc_length(float *x2d, float *z2d, int nx,
                    int nz, float *arc1_len, float *arc2_len)
{
  size_t iptr1, iptr2, iptr3, iptr4; 
  float x_len, z_len, dh_len;
  for(int k=1; k<nz; k++)
  {
    iptr1 = k*nx + 0;
    iptr2 = (k-1)*nx + 0;
    x_len = x2d[iptr1] - x2d[iptr2];
    z_len = z2d[iptr1] - z2d[iptr2];
    dh_len = sqrt(pow(x_len,2) + pow(z_len,2));
    arc1_len[k] = arc1_len[k-1] + dh_len;

    iptr3 = k*nx + nx-1;
    iptr4 = (k-1)*nx + nx-1;
    x_len = x2d[iptr3] - x2d[iptr4];
    z_len = z2d[iptr3] - z2d[iptr4];
    dh_len = sqrt(pow(x_len,2) + pow(z_len,2));
    arc2_len[k] = arc2_len[k-1] + dh_len;
  }

  return 0;
}
