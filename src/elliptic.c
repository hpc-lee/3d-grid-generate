#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "elliptic.h"
#include "constants.h"
#include "solver.h"
#include "lib_mem.h"

int
diri_gene(gd_t *gdcurv, par_t *par)
{
  float err = par->i_err;
  int max_iter = par->max_iter;
  float coef = par->coef;

  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float *P = NULL; //source term P
  float *Q = NULL; //source term Q
  P = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source P");
  Q = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source Q");

  float *x2d_tmp =NULL;
  float *z2d_tmp =NULL;
  x2d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "x2d_tmp");
  z2d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "z2d_tmp");

  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  size_t iptr, iptr1, iptr2;
  float resi, resk;
  float max_resi, max_resk;
  
  // before update grid. use init grid to
  // calculate ghost points and bdry g11,g22 
  float *p_x; // point_x_bdry
  float *p_z; // point_z_bdry
  float *g11_x;  // g11_x_bdry
  float *g22_z;  // g22_z_bdry
  p_x = (float *)mem_calloc_1d_float(nz*2*2, 0.0, "p_x");
  p_z = (float *)mem_calloc_1d_float(nx*2*2, 0.0, "p_z");
  g11_x = (float *)mem_calloc_1d_float(nz*2, 0.0, "g11_x"); 
  g22_z = (float *)mem_calloc_1d_float(nx*2, 0.0, "g22_z"); 
  ghost_cal(x2d,z2d,nx,nz,p_x,p_z,g11_x,g22_z);

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr = k*nx + i;
      x2d_tmp[iptr] = x2d[iptr];
      z2d_tmp[iptr] = z2d[iptr];
    }
  }

  int flag_true = 1;
  while(flag_true)
  {
    // update solver
    update_SOR(x2d,z2d,x2d_tmp,z2d_tmp,nx,nz,P,Q,w);
    Niter += 1;

    max_resi = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*nx + i;
        iptr1 = k*nx + i+1;
        iptr2 = (k+1)*nx + i;
        resi = fabs((x2d_tmp[iptr] - x2d[iptr])/(x2d[iptr1]-x2d[iptr]));
        resk = fabs((z2d_tmp[iptr] - z2d[iptr])/(z2d[iptr2]-z2d[iptr]));
        max_resi = fmax(max_resi,resi);
        max_resk = fmax(max_resk,resk);
      }
    }

    // copy coord
    for(int k=0; k<nz; k++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*nx + i;
        x2d[iptr] = x2d_tmp[iptr];
        z2d[iptr] = z2d_tmp[iptr];
      }
    }
    
    set_src_diri(x2d,z2d,nx,nz,P,Q,p_x,p_z,g11_x,g22_z,coef);

    fprintf(stdout,"number of iter is %d\n", Niter);
    fprintf(stdout,"max_resi is %f, max_resk is %f\n", max_resi, max_resk);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resk < err) {
      flag_true = 0;
    }

  }

  free(x2d_tmp);
  free(z2d_tmp);
  free(P);
  free(Q);
  free(p_x);
  free(p_z);
  free(g11_x);
  free(g22_z);

  return 0;
}

int
set_src_diri(float *x2d, float *z2d, int nx, int nz, 
             float *P, float *Q, float *p_x, float *p_z,
             float *g11_x, float *g22_z, float coef)
{
  float *p_x1;
  float *p_x2;
  float *p_z1;
  float *p_z2;
  p_x1 = p_x;
  p_x2 = p_x + nz*2;
  p_z1 = p_z;
  p_z2 = p_z + nx*2;

  size_t iptr, iptr1, iptr2, iptr3;
  float x_xi,z_xi,x_xixi,z_xixi;
  float x_zt,z_zt,x_ztzt,z_ztzt;
  float g11, g22;

  // bdry z1 zt=0
  for(int i=1; i<nx-1; i++)
  {
    iptr  = i;         //(i,0)
    iptr1 = i+1;       //(i+1,0)
    iptr2 = i-1;       //(i-1,0)
    x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
    z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]);
    x_xixi = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
    z_xixi = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

    iptr3 = 1*nx + i;  //(i,1)
    x_zt = x2d[iptr3] - x2d[iptr];
    z_zt = z2d[iptr3] - z2d[iptr];
    x_ztzt = p_z1[i] + x2d[iptr3] - 2*x2d[iptr];
    z_ztzt = p_z1[i+nx] + z2d[iptr3] - 2*z2d[iptr];
    g11 = pow(x_xi,2) + pow(z_xi,2);

    P[iptr] = -(x_xi*x_xixi + z_xi*z_xixi)/g11  
              -(x_xi*x_ztzt + z_xi*z_ztzt)/g22_z[i];
    Q[iptr] = -(x_zt*x_xixi + z_zt*z_xixi)/g11  
              -(x_zt*x_ztzt + z_zt*z_ztzt)/g22_z[i];
  }

  // bdry z2 zt=1
  for(int i=1; i<nx-1; i++)
  {
    iptr  = (nz-1)*nx + i;      //(i,nz-1)
    iptr1 = (nz-1)*nx + (i+1);  //(i+1,nz-1)
    iptr2 = (nz-1)*nx + (i-1);  //(i-1,nz-1)
    x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
    z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]);
    x_xixi = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
    z_xixi = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

    iptr3 = (nz-2)*nx + i;  //(i,nz-2)
    x_zt = x2d[iptr] - x2d[iptr3];
    z_zt = z2d[iptr] - z2d[iptr3];
    x_ztzt = p_z2[i] + x2d[iptr3] - 2*x2d[iptr];
    z_ztzt = p_z2[i+nx] + z2d[iptr3] - 2*z2d[iptr];
    g11 = pow(x_xi,2) + pow(z_xi,2);

    P[iptr] = -(x_xi*x_xixi + z_xi*z_xixi)/g11  
              -(x_xi*x_ztzt + z_xi*z_ztzt)/g22_z[i+nx];
    Q[iptr] = -(x_zt*x_xixi + z_zt*z_xixi)/g11  
              -(x_zt*x_ztzt + z_zt*z_ztzt)/g22_z[i+nx];
  }
 
  // bdry x1 xi=0
  for(int k=1; k<nz-1; k++)
  {
    iptr  = k*nx;         //(0,k)
    iptr1 = (k+1)*nx;     //(0,k+1)
    iptr2 = (k-1)*nx;     //(0,k-1)
    x_zt = 0.5*(x2d[iptr1] - x2d[iptr2]);
    z_zt = 0.5*(z2d[iptr1] - z2d[iptr2]);
    x_ztzt = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
    z_ztzt = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

    iptr3 = k*nx + 1;  //(1,k)
    x_xi = x2d[iptr3] - x2d[iptr];
    z_xi = z2d[iptr3] - z2d[iptr];
    x_xixi = p_x1[k] + x2d[iptr3] - 2*x2d[iptr];
    z_xixi = p_x1[k+nz] + z2d[iptr3] - 2*z2d[iptr];
    g22 = pow(x_zt,2) + pow(z_zt,2);

    P[iptr] = -(x_xi*x_xixi + z_xi*z_xixi)/g11_x[k]  
              -(x_xi*x_ztzt + z_xi*z_ztzt)/g22;
    Q[iptr] = -(x_zt*x_xixi + z_zt*z_xixi)/g11_x[k]  
              -(x_zt*x_ztzt + z_zt*z_ztzt)/g22;
  }

  // bdry x2 xi=1
  for(int k=1; k<nz-1; k++)
  {
    iptr  = k*nx + (nx-1);         //(nx-1,k)
    iptr1 = (k+1)*nx + (nx-1);     //(nx-1,k+1)
    iptr2 = (k-1)*nx + (nx-1);     //(nx-1,k-1)
    x_zt = 0.5*(x2d[iptr1] - x2d[iptr2]);
    z_zt = 0.5*(z2d[iptr1] - z2d[iptr2]);
    x_ztzt = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
    z_ztzt = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

    iptr3 = k*nx + (nx-2);  //(nx-2,k)
    x_xi = x2d[iptr] - x2d[iptr3];
    z_xi = z2d[iptr] - z2d[iptr3];
    x_xixi = p_x2[k] + x2d[iptr3] - 2*x2d[iptr];
    z_xixi = p_x2[k+nz] + z2d[iptr3] - 2*z2d[iptr];
    g22 = pow(x_zt,2) + pow(z_zt,2);

    P[iptr] = -(x_xi*x_xixi + z_xi*z_xixi)/g11_x[k+nz]  
              -(x_xi*x_ztzt + z_xi*z_ztzt)/g22;
    Q[iptr] = -(x_zt*x_xixi + z_zt*z_xixi)/g11_x[k+nz]  
              -(x_zt*x_ztzt + z_zt*z_ztzt)/g22;
  }
  // first use bdry z1 z2 to interp inner point
  float xi,zt,c0,c1,r0,r1;
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      zt = (1.0*k)/(nz-1);
      c0 = 1-zt;
      c1 = zt;

      r0 = exp(coef*zt);
      r1 = exp(coef*(1-zt)); 
      
      iptr  = k*nx + i;
      iptr1 = i;
      iptr2 = (nz-1)*nx + i;
      P[iptr] = r0*P[iptr1] + r1*P[iptr2];
      Q[iptr] = c0*Q[iptr1] + c1*Q[iptr2];
    }
  }
  
  // then use bdry x1 x2 to interp inner point
  // NOTE: point (:,0:4),(:,nz-5:nz-1) P Q calculate
  // only by z1,z2
  // point(:,5:nz-6), P Q calculte by z1 z2 x1 x2
  // so need superposition
  for(int k=5; k<nz-5; k++) {
    for(int i=1; i<nx-1; i++)
    {
      xi = (1.0*i)/(nx-1);

      c0 = 1-xi;
      c1 = xi;

      r0 = exp(coef*xi);
      r1 = exp(coef*(1-xi)); 
      
      iptr  = k*nx + i;
      iptr1 = k*nx + 0;
      iptr2 = k*nx + (nx-1);
      P[iptr] = P[iptr] + c0*P[iptr1] + c1*P[iptr2];
      //P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
      Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
    }
  }

  return 0;
}

int
ghost_cal(float *x2d, float *z2d, int nx, int nz, float *p_x, float *p_z,
          float *g11_x, float *g22_z)
{
  
  size_t iptr1,iptr2,iptr3,iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float x_xi0,z_xi0,x_zt0,z_zt0;
  float vn_xi, vn_zt, vn_xi0, vn_zt0;
  float len_vn, dot;
  float *p_x1;
  float *p_x2;
  float *p_z1;
  float *p_z2;
  p_x1 = p_x;
  p_x2 = p_x + nz*2;
  p_z1 = p_z;
  p_z2 = p_z + nx*2;
  // bdry x1
  for(int k=1; k<nz-1; k++)
  {
    iptr1 = k*nx + 1;    //(1,k)
    iptr2 = k*nx + 0;    //(0,k)
    iptr3 = (k+1)*nx + 0;//(0,k+1)
    iptr4 = (k-1)*nx + 0;//(0,k-1)
    x_xi0 = x2d[iptr1] - x2d[iptr2];
    z_xi0 = z2d[iptr1] - z2d[iptr2];
    x_zt = 0.5*(x2d[iptr3] - x2d[iptr4]);
    z_zt = 0.5*(z2d[iptr3] - z2d[iptr4]);
    // orth vector
    vn_xi =  z_zt;
    vn_zt = -x_zt;
    len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
    // norm
    vn_xi0 = vn_xi/len_vn;
    vn_zt0 = vn_zt/len_vn;
    // projection from r_xi0 to vn
    dot = x_xi0*vn_xi0 + z_xi0*vn_zt0;
    x_xi = dot*vn_xi0;
    z_xi = dot*vn_zt0;

    p_x1[k] = x2d[iptr2] - x_xi;
    p_x1[k+nz] = z2d[iptr2] - z_xi;
    g11_x[k] = pow(x_xi,2) + pow(z_xi,2);
  }
  // bdry x2
  for(int k=1; k<nz-1; k++)
  {
    iptr1 = k*nx + (nx-1);    //(nx-1,k)
    iptr2 = k*nx + (nx-2);    //(nx-2,k)
    iptr3 = (k+1)*nx + (nx-1);//(nx-1,k+1)
    iptr4 = (k-1)*nx + (nx-1);//(nx-1,k-1)
    x_xi0 = x2d[iptr1] - x2d[iptr2];
    z_xi0 = z2d[iptr1] - z2d[iptr2];
    x_zt = 0.5*(x2d[iptr3] - x2d[iptr4]);
    z_zt = 0.5*(z2d[iptr3] - z2d[iptr4]);
    // orth vector
    vn_xi =  z_zt;
    vn_zt = -x_zt;
    len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
    // norm
    vn_xi0 = vn_xi/len_vn;
    vn_zt0 = vn_zt/len_vn;
    // projection from r_xi0 to vn
    dot = x_xi0*vn_xi0 + z_xi0*vn_zt0;
    x_xi = dot*vn_xi0;
    z_xi = dot*vn_zt0;

    p_x2[k] = x2d[iptr1] + x_xi;
    p_x2[k+nz] = z2d[iptr1] + z_xi;
    g11_x[k+nz] = pow(x_xi,2) + pow(z_xi,2);
  }
  // bdry z1
  for(int i=1; i<nx-1; i++)
  {
    iptr1 = 1*nx + i;    //(i,1)
    iptr2 = 0*nx + i;    //(i,0)
    iptr3 = 0*nx + i+1;  //(i+1,0)
    iptr4 = 0*nx + i-1;  //(i-1,0)
    x_zt0 = x2d[iptr1] - x2d[iptr2];
    z_zt0 = z2d[iptr1] - z2d[iptr2];
    x_xi = 0.5*(x2d[iptr3] - x2d[iptr4]);
    z_xi = 0.5*(z2d[iptr3] - z2d[iptr4]);
    // orth vector
    vn_xi = -z_xi;
    vn_zt =  x_xi;
    len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
    // norm
    vn_xi0 = vn_xi/len_vn;
    vn_zt0 = vn_zt/len_vn;
    // projection from r_zt0 to vn
    dot = x_zt0*vn_xi0 + z_zt0*vn_zt0;
    x_zt = dot*vn_xi0;
    z_zt = dot*vn_zt0;

    p_z1[i] = x2d[iptr2] - x_zt;
    p_z1[i+nx] = z2d[iptr2] - z_zt;
    g22_z[i] = pow(x_zt,2) + pow(z_zt,2);
  }
  // bdry z2
  for(int i=1; i<nx-1; i++)
  {
    iptr1 = (nz-1)*nx + i;    //(i,nz-1)
    iptr2 = (nz-2)*nx + i;    //(i,nz-2)
    iptr3 = (nz-1)*nx + i+1;  //(i+1,nz-1)
    iptr4 = (nz-1)*nx + i-1;  //(i-1,nz-1)
    x_zt0 = x2d[iptr1] - x2d[iptr2];
    z_zt0 = z2d[iptr1] - z2d[iptr2];
    x_xi = 0.5*(x2d[iptr3] - x2d[iptr4]);
    z_xi = 0.5*(z2d[iptr3] - z2d[iptr4]);
    // orth vector
    vn_xi = -z_xi;
    vn_zt =  x_xi;
    len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
    // norm
    vn_xi0 = vn_xi/len_vn;
    vn_zt0 = vn_zt/len_vn;
    // projection from r_zt0 to vn
    dot = x_zt0*vn_xi0 + z_zt0*vn_zt0;
    x_zt = dot*vn_xi0;
    z_zt = dot*vn_zt0;

    p_z2[i] = x2d[iptr1] + x_zt;
    p_z2[i+nx] = z2d[iptr1] + z_zt;
    g22_z[i+nx] = pow(x_zt,2) + pow(z_zt,2);
  }

  return 0;
}

int
higen_gene(gd_t *gdcurv, par_t *par)
{
  float err = par->i_err;
  int max_iter = par->max_iter;
  float coef = par->coef;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float *P = NULL; //source term P
  float *Q = NULL; //source term Q
  P = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source P");
  Q = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source Q");

  float *x2d_tmp =NULL;
  float *z2d_tmp =NULL;
  x2d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "x2d_tmp");
  z2d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "z2d_tmp");

  float dx1,dx2,dz1,dz2;
  if(par->dire_itype == Z_DIRE)
  {
    dx1 = par->distance[0];
    dx2 = par->distance[1];
    dz1 = par->distance[2];
    dz2 = par->distance[3];
  }
  if(par->dire_itype == X_DIRE)
  {
    dz1 = par->distance[0];
    dz2 = par->distance[1];
    dx1 = par->distance[2];
    dx2 = par->distance[3];
  }
  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  size_t iptr, iptr1, iptr2;
  float resi, resk;
  float max_resi, max_resk;

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr = k*nx + i;
      x2d_tmp[iptr] = x2d[iptr];
      z2d_tmp[iptr] = z2d[iptr];
    }
  }

  int flag_true = 1;
  while(flag_true)
  {
    // update solver
    update_SOR(x2d,z2d,x2d_tmp,z2d_tmp,nx,nz,P,Q,w);
    Niter += 1;

    max_resi = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr = k*nx + i;
        iptr1 = k*nx + i+1;
        iptr2 = (k+1)*nx + i;
        resi = fabs((x2d_tmp[iptr] - x2d[iptr])/(x2d[iptr1]-x2d[iptr]));
        resk = fabs((z2d_tmp[iptr] - z2d[iptr])/(z2d[iptr2]-z2d[iptr]));
        max_resi = fmax(max_resi,resi);
        max_resk = fmax(max_resk,resk);
      }
    }

    // copy coord
    for(int k=0; k<nz; k++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*nx + i;
        x2d[iptr] = x2d_tmp[iptr];
        z2d[iptr] = z2d_tmp[iptr];
      }
    }
    
    set_src_higen(x2d,z2d,nx,nz,P,Q,coef,dx1,dx2,dz1,dz2);

    fprintf(stdout,"number of iter is %d\n", Niter);
    fprintf(stdout,"max_resi is %f, max_resk is %f\n", max_resi, max_resk);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resk < err) {
      flag_true = 0;
    }

  }

  free(x2d_tmp);
  free(z2d_tmp);
  free(P);
  free(Q);

  return 0;
}

int
set_src_higen(float *x2d, float *z2d, int nx, int nz,
              float *P, float *Q, float coef,
              float dx1, float dx2, float dz1, float dz2)
{
  float theta0 = PI/2;
  float a = 0.1;
  size_t iptr,iptr1,iptr2,iptr3;
  float dot, len_xi, len_zt, dif_dis;
  float cos_theta, theta, dif_theta;
  float x_xi, z_xi, x_zt, z_zt;
  // bdry z1 zt=0
  for(int i=1; i<nx-1; i++)
  {
    iptr  = i;         //(i,0)
    iptr1 = i+1;       //(i+1,0)
    iptr2 = i-1;       //(i-1,0)
    iptr3 = 1*nx + i;  //(i,1)

    x_xi = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_xi = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_zt = x2d[iptr3] - x2d[iptr];
    z_zt = z2d[iptr3] - z2d[iptr];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    P[iptr] = P[iptr] - a*tanh(dif_theta);

    dif_dis = (dz1-len_zt)/dz1;
    Q[iptr] = Q[iptr] + a*tanh(dif_dis);
  }

  // NOTE: P and Q update need change sign
  // P is +a and Q is -a
  // bdry z2 zt=0
  for(int i=1; i<nx-1; i++)
  {
    iptr  = (nz-1)*nx + i;      //(i,nz-1)
    iptr1 = (nz-1)*nx + i+1;    //(i+1,nz-1)
    iptr2 = (nz-1)*nx + i-1;    //(i-1,nz-1)
    iptr3 = (nz-2)*nx + i;      //(i,nz-2)

    x_xi = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_xi = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_zt = x2d[iptr] - x2d[iptr3];
    z_zt = z2d[iptr] - z2d[iptr3];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    P[iptr] = P[iptr] + a*tanh(dif_theta);

    dif_dis = (dz2-len_zt)/dz2;
    Q[iptr] = Q[iptr] - a*tanh(dif_dis);
  }

  // bdry x1 xi=0
  for(int k=1; k<nz-1; k++)
  {
    iptr  =  k*nx + 0;         //(0,k)
    iptr1 = (k+1)*nx + 0;      //(0,k+1)
    iptr2 = (k-1)*nx + 0;      //(0,k-1)
    iptr3 = k*nx + 1;  //(1,k)

    x_zt = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_zt = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_xi = x2d[iptr3] - x2d[iptr];
    z_xi = z2d[iptr3] - z2d[iptr];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    Q[iptr] = Q[iptr] - a*tanh(dif_theta);

    dif_dis = (dx1-len_xi)/dx1;
    P[iptr] = P[iptr] + a*tanh(dif_dis);
  }
  // bdry x2 xi=1
  for(int k=1; k<nz-1; k++)
  {
    iptr  =  k*nx + (nx-1);         //(nx-1,k)
    iptr1 = (k+1)*nx + (nx-1);      //(nx-1,k+1)
    iptr2 = (k-1)*nx + (nx-1);      //(nx-1,k-1)
    iptr3 = k*nx + (nx-2);  //(nx-2,k)

    x_zt = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_zt = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_xi = x2d[iptr] - x2d[iptr3];
    z_xi = z2d[iptr] - z2d[iptr3];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    Q[iptr] = Q[iptr] + a*tanh(dif_theta);

    dif_dis = (dx2-len_xi)/dx2;
    P[iptr] = P[iptr] - a*tanh(dif_dis);
  }

  // first use bdry z1 z2 to interp inner point
  float xi,zt,c0,c1,r0,r1;
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      zt = (1.0*k)/(nz-1);
      c0 = 1-zt;
      c1 = zt;

      r0 = exp(coef*zt);
      r1 = exp(coef*(1-zt)); 
      
      iptr  = k*nx + i;
      iptr1 = i;
      iptr2 = (nz-1)*nx + i;
      P[iptr] = r0*P[iptr1] + r1*P[iptr2];
      Q[iptr] = c0*Q[iptr1] + c1*Q[iptr2];
    }
  }
  // then use bdry x1 x2 to interp inner point
  // NOTE: point (:,0:4),(:,nz-5:nz-1) P Q calculate
  // only by z1,z2
  // point(:,5:nz-6), P Q calculte by z1 z2 x1 x2
  // so need superposition
  for(int k=5; k<nz-5; k++) {
    for(int i=1; i<nx-1; i++)
    {
      xi = (1.0*i)/(nx-1);

      c0 = 1-xi;
      c1 = xi;

      r0 = exp(coef*xi);
      r1 = exp(coef*(1-xi)); 
      
      iptr  = k*nx + i;
      iptr1 = k*nx + 0;
      iptr2 = k*nx + (nx-1);
      P[iptr] = P[iptr] + c0*P[iptr1] + c1*P[iptr2];
      //P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
      Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
    }
  }

  return 0;
}
