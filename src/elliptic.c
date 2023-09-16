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
  float err = par->iter_err;
  int max_iter = par->max_iter;
  float coef = par->coef;
  int first_dire_itype = par->first_dire_itype;
  int second_dire_itype = par->second_dire_itype;

  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  float *P = NULL; //source term P
  float *Q = NULL; //source term Q
  float *R = NULL; //source term R
  P = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source P");
  Q = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source Q");
  R = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source R");

  float *x3d_tmp =NULL;
  float *y3d_tmp =NULL;
  float *z3d_tmp =NULL;
  x3d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "x3d_tmp");
  y3d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "y3d_tmp");
  z3d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "z3d_tmp");

  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  size_t iptr, iptr1, iptr2, iptr3;
  float resi, resj, resk;
  float max_resi, max_resj, max_resk;
  
  // before update grid. use init grid to
  // calculate ghost points and bdry g11,g22,g33 
  float *p_x; // point_x_bdry
  float *p_y; // point_y_bdry
  float *p_z; // point_z_bdry
  float *g11_x;  // g11_x_bdry
  float *g22_y;  // g22_y_bdry
  float *g33_z;  // g33_z_bdry
  p_x = (float *)mem_calloc_1d_float(nz*ny*3*2, 0.0, "p_x");
  p_y = (float *)mem_calloc_1d_float(nz*nx*3*2, 0.0, "p_y");
  p_z = (float *)mem_calloc_1d_float(ny*nx*3*2, 0.0, "p_z");
  g11_x = (float *)mem_calloc_1d_float(nz*ny*2, 0.0, "g11_x");
  g22_y = (float *)mem_calloc_1d_float(nz*nx*2, 0.0, "g22_y");
  g33_z = (float *)mem_calloc_1d_float(ny*nx*2, 0.0, "g33_z");

  ghost_cal(x3d,y3d,z3d,nx,ny,nz,p_x,p_y,p_z,g11_x,g22_y,g33_z);

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + j*siz_iy + i;
        x3d_tmp[iptr] = x3d[iptr];
        y3d_tmp[iptr] = y3d[iptr];
        z3d_tmp[iptr] = z3d[iptr];
      }
    }
  }

  int flag_true = 1;
  while(flag_true)
  {
    // update solver
    update_SOR(x3d,y3d,z3d,x3d_tmp,y3d_tmp,z3d_tmp,nx,ny,nz,P,Q,R,w);
    Niter += 1;

    max_resi = 0.0;
    max_resj = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + i+1;
          iptr2 = k*siz_iz + (j+1)*siz_iy + i;
          iptr3 = (k+1)*siz_iz + j*siz_iy + i;
          resi = fabs((x3d_tmp[iptr] - x3d[iptr])/(x3d[iptr1]-x3d[iptr]));
          resj = fabs((y3d_tmp[iptr] - y3d[iptr])/(y3d[iptr2]-y3d[iptr]));
          resk = fabs((z3d_tmp[iptr] - z3d[iptr])/(z3d[iptr3]-z3d[iptr]));
          max_resi = fmax(max_resi,resi);
          max_resj = fmax(max_resj,resj);
          max_resk = fmax(max_resk,resk);
        }
      }
    }

    // copy coord
    for(int k=0; k<nz; k++) {
      for(int j=0; j<ny; j++) {
        for(int i=0; i<nx; i++)
        {
          iptr = k*siz_iz + j*siz_iy + i;
          x3d[iptr] = x3d_tmp[iptr];
          y3d[iptr] = y3d_tmp[iptr];
          z3d[iptr] = z3d_tmp[iptr];
        }
      }
    }
    
    set_src_diri(x3d,y3d,z3d,nx,ny,nz,P,Q,R,p_x,p_y,p_z,g11_x,g22_y,g33_z,coef,
                 first_dire_itype,second_dire_itype);

    fprintf(stdout,"number of iter is %d\n", Niter);
    fprintf(stdout,"max_resi is %f, max_resj is %f, max_resk is %f\n", max_resi, max_resj, max_resk);
    fflush(stdout);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resj < err && max_resk < err) {
      flag_true = 0;
    }

  }

  free(x3d_tmp);
  free(y3d_tmp);
  free(z3d_tmp);
  free(P);
  free(Q);
  free(R);
  free(p_x);
  free(p_y);
  free(p_z);
  free(g11_x);
  free(g22_y);
  free(g33_z);

  return 0;
}

int
set_src_diri(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, 
             float *P, float *Q, float *R, float *p_x, float *p_y, 
             float *p_z, float *g11_x, float *g22_y, float *g33_z, 
             float coef, int first_dire_itype, int second_dire_itype)
{
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float temp_x,temp_y,temp_z; 
  float temp1,temp2,temp3; 
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xixi,y_xixi,z_xixi;
  float x_etet,y_etet,z_etet;
  float x_ztzt,y_ztzt,z_ztzt;
  float x_xizt,y_xizt,z_xizt;
  float x_xiet,y_xiet,z_xiet;
  float x_etzt,y_etzt,z_etzt;
  float g11,g22,g33,g12,g13,g23;
  float c1,c2,c3,c4;
  float *p_x1,*p_x2;
  float *p_y1,*p_y2;
  float *p_z1,*p_z2;
  float *g11_x1, *g11_x2;
  float *g22_y1, *g22_y2;
  float *g33_z1, *g33_z2;
  size_t size;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;
  p_x1 = p_x;
  p_x2 = p_x + ny*nz*3;
  p_y1 = p_y;
  p_y2 = p_y + nx*nz*3;
  p_z1 = p_z;
  p_z2 = p_z + nx*ny*3;
  g11_x1 = g11_x;
  g11_x2 = g11_x + ny*nz;
  g22_y1 = g22_y;
  g22_y2 = g22_y + nx*nz;
  g33_z1 = g33_z;
  g33_z2 = g33_z + nx*ny;

  // bdry x1 
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr  = k*siz_iz + j*siz_iy + 0;     //(0,j,k)
      iptr1 = k*siz_iz + (j+1)*siz_iy + 0; //(0,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + 0; //(0,j-1,k)
      iptr3 = (k+1)*siz_iz + j*siz_iy + 0; //(0,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + 0; //(0,j,k-1)
      x_et = 0.5*(x3d[iptr1] - x3d[iptr2]); 
      y_et = 0.5*(y3d[iptr1] - y3d[iptr2]); 
      z_et = 0.5*(z3d[iptr1] - z3d[iptr2]); 
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]); 
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]); 
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]); 
      x_etet = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
      y_etet = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
      z_etet = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
      x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
      y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
      z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

      iptr1 = (k-1)*siz_iz + (j-1)*siz_iy + 0; //(0,j-1,k-1)
      iptr2 = (k-1)*siz_iz + (j+1)*siz_iy + 0; //(0,j+1,k-1)
      iptr3 = (k+1)*siz_iz + (j+1)*siz_iy + 0; //(0,j+1,k+1)
      iptr4 = (k+1)*siz_iz + (j-1)*siz_iy + 0; //(0,j-1,k+1)

      x_etzt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
      y_etzt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
      z_etzt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
      
      iptr1 = k*siz_iz + j*siz_iy + 1;     //(1,j,k)
      iptr2 = k*ny + j;
      size = ny*nz;
      x_xi = x3d[iptr1] - x3d[iptr];
      y_xi = y3d[iptr1] - y3d[iptr];
      z_xi = z3d[iptr1] - z3d[iptr];
      x_xixi = p_x1[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
      y_xixi = p_x1[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
      z_xixi = p_x1[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

      g22 = x_et*x_et + y_et*y_et + z_et*z_et;
      g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
      g23 = x_et*x_zt + y_et*y_zt + z_et*z_zt;

      c1 = 1/(g22*g33 - g23*g23); // 1/alpha1
      c2 = g23/g33;
      c3 = g23/g22;
      c4 = 1/g11_x1[iptr2];

      temp1 = -c4*(x_xi*x_xixi + y_xi*y_xixi + z_xi*z_xixi);
      temp2 = -c4*((x_et-c2*x_zt)*x_xixi + (y_et-c2*y_zt)*y_xixi + (z_et-c2*z_zt)*z_xixi);
      temp3 = -c4*((x_zt-c3*x_et)*x_xixi + (y_zt-c3*y_et)*y_xixi + (z_zt-c3*z_et)*z_xixi);
      temp_x = g33*x_etet + g22*x_ztzt - 2*g23*x_etzt;
      temp_y = g33*y_etet + g22*y_ztzt - 2*g23*y_etzt;
      temp_z = g33*z_etet + g22*z_ztzt - 2*g23*z_etzt;

      P[iptr] = temp1 - c1*(x_xi*temp_x + y_xi*temp_y + z_xi*temp_z);

      Q[iptr] = temp2 - c1*((x_et-c2*x_zt)*temp_x + (y_et-c2*y_zt)*temp_y 
              + (z_et-c2*z_zt)*temp_z);

      R[iptr] = temp3 - c1*((x_zt-c3*x_et)*temp_x + (y_zt-c3*y_et)*temp_y 
              + (z_zt-c3*z_et)*temp_z);
    }
  }

  // bdry x2 
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr  = k*siz_iz + j*siz_iy + nx-1;     //(nx-1,j,k)
      iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1; //(nx-1,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1; //(nx-1,j-1,k)
      iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k-1)
      x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
      x_etet = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
      y_etet = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
      z_etet = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
      x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
      y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
      z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

      iptr1 = (k-1)*siz_iz + (j-1)*siz_iy + nx-1; //(nx-1,j-1,k-1)
      iptr2 = (k-1)*siz_iz + (j+1)*siz_iy + nx-1; //(nx-1,j+1,k-1)
      iptr3 = (k+1)*siz_iz + (j+1)*siz_iy + nx-1; //(nx-1,j+1,k+1)
      iptr4 = (k+1)*siz_iz + (j-1)*siz_iy + nx-1; //(nx-1,j-1,k+1)

      x_etzt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
      y_etzt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
      z_etzt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
      
      iptr1 = k*siz_iz + j*siz_iy + nx-2;     //(nx-2,j,k)
      iptr2 = k*ny + j;
      size = ny*nz;
      x_xi = x3d[iptr] - x3d[iptr1];
      y_xi = y3d[iptr] - y3d[iptr1];
      z_xi = z3d[iptr] - z3d[iptr1];
      x_xixi = p_x2[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
      y_xixi = p_x2[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
      z_xixi = p_x2[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

      g22 = x_et*x_et + y_et*y_et + z_et*z_et;
      g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
      g23 = x_et*x_zt + y_et*y_zt + z_et*z_zt;

      c1 = 1/(g22*g33 - g23*g23); // 1/alpha1
      c2 = g23/g33;
      c3 = g23/g22;
      c4 = 1/g11_x2[iptr2];

      temp1 = -c4*(x_xi*x_xixi + y_xi*y_xixi + z_xi*z_xixi);
      temp2 = -c4*((x_et-c2*x_zt)*x_xixi + (y_et-c2*y_zt)*y_xixi + (z_et-c2*z_zt)*z_xixi);
      temp3 = -c4*((x_zt-c3*x_et)*x_xixi + (y_zt-c3*y_et)*y_xixi + (z_zt-c3*z_et)*z_xixi);
      temp_x = g33*x_etet + g22*x_ztzt - 2*g23*x_etzt;
      temp_y = g33*y_etet + g22*y_ztzt - 2*g23*y_etzt;
      temp_z = g33*z_etet + g22*z_ztzt - 2*g23*z_etzt;

      P[iptr] = temp1 - c1*(x_xi*temp_x + y_xi*temp_y + z_xi*temp_z);

      Q[iptr] = temp2 - c1*((x_et-c2*x_zt)*temp_x + (y_et-c2*y_zt)*temp_y 
              + (z_et-c2*z_zt)*temp_z);

      R[iptr] = temp3 - c1*((x_zt-c3*x_et)*temp_x + (y_zt-c3*y_et)*temp_y 
              + (z_zt-c3*z_et)*temp_z);
    }
  }
  // bdry y1 
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = k*siz_iz + 0*siz_iy + i;     //(i,0,k)
      iptr1 = k*siz_iz + 0*siz_iy + i+1;   //(i+1,0,k)
      iptr2 = k*siz_iz + 0*siz_iy + i-1;   //(i-1,0,k)
      iptr3 = (k+1)*siz_iz + 0*siz_iy + i; //(i,0,k+1)
      iptr4 = (k-1)*siz_iz + 0*siz_iy + i; //(i,0,k-1)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]); 
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]); 
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]); 
      x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
      y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
      z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
      x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
      y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
      z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

      iptr1 = (k-1)*siz_iz + 0*siz_iy + i-1; //(i-1,0,k-1)
      iptr2 = (k-1)*siz_iz + 0*siz_iy + i+1; //(i+1,0,k-1)
      iptr3 = (k+1)*siz_iz + 0*siz_iy + i+1; //(i+1,0,k+1)
      iptr4 = (k+1)*siz_iz + 0*siz_iy + i-1; //(i-1,0,k+1)

      x_xizt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
      y_xizt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
      z_xizt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
      
      iptr1 = k*siz_iz + 1*siz_iy + i;     //(i,1,k)
      iptr2 = k*nx + i;
      size = nx*nz;
      x_et = x3d[iptr1] - x3d[iptr];
      y_et = y3d[iptr1] - y3d[iptr];
      z_et = z3d[iptr1] - z3d[iptr];
      x_etet = p_y1[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
      y_etet = p_y1[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
      z_etet = p_y1[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

      g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
      g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
      g13 = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;

      c1 = 1/(g11*g33 - g13*g13); // 1/alpha2
      c2 = g13/g33;
      c3 = g13/g11;
      c4 = 1/g22_y1[iptr2];

      temp1 = -c4*((x_xi-c2*x_zt)*x_etet + (y_xi-c2*y_zt)*y_etet + (z_xi-c2*z_zt)*z_etet);
      temp2 = -c4*(x_et*x_etet + y_et*y_etet + z_et*z_etet);
      temp3 = -c4*((x_zt-c3*x_xi)*x_etet + (y_zt-c3*y_xi)*y_etet + (z_zt-c3*z_xi)*z_etet);
      temp_x = g33*x_xixi + g11*x_ztzt - 2*g13*x_xizt;
      temp_y = g33*y_xixi + g11*y_ztzt - 2*g13*y_xizt;
      temp_z = g33*z_xixi + g11*z_ztzt - 2*g13*z_xizt;

      P[iptr] = temp1 - c1*((x_xi-c2*x_zt)*temp_x + (y_xi-c2*y_zt)*temp_y
                          + (z_xi-c2*z_zt)*temp_z);

      Q[iptr] = temp2 - c1*(x_et*temp_x + y_et*temp_y + z_et*temp_z);

      R[iptr] = temp3 - c1*((x_zt-c3*x_xi)*temp_x + (y_zt-c3*y_xi)*temp_y 
                          + (z_zt-c3*z_xi)*temp_z);
    }
  }
  // bdry y2 
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = k*siz_iz + (ny-1)*siz_iy + i;     //(i,ny-1,k)
      iptr1 = k*siz_iz + (ny-1)*siz_iy + i+1;   //(i+1,ny-1,k)
      iptr2 = k*siz_iz + (ny-1)*siz_iy + i-1;   //(i-1,ny-1,k)
      iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k+1)
      iptr4 = (k-1)*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k-1)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]); 
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]); 
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]); 
      x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
      y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
      z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
      x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
      y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
      z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

      iptr1 = (k-1)*siz_iz + (ny-1)*siz_iy + i-1; //(i-1,ny-1,k-1)
      iptr2 = (k-1)*siz_iz + (ny-1)*siz_iy + i+1; //(i+1,ny-1,k-1)
      iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i+1; //(i+1,ny-1,k+1)
      iptr4 = (k+1)*siz_iz + (ny-1)*siz_iy + i-1; //(i-1,ny-1,k+1)

      x_xizt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
      y_xizt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
      z_xizt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
      
      iptr1 = k*siz_iz + (ny-2)*siz_iy + i;     //(i,1,k)
      iptr2 = k*nx + i;
      size = nx*nz;
      x_et = x3d[iptr] - x3d[iptr1];
      y_et = y3d[iptr] - y3d[iptr1];
      z_et = z3d[iptr] - z3d[iptr1];
      x_etet = p_y2[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
      y_etet = p_y2[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
      z_etet = p_y2[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

      g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
      g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
      g13 = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;

      c1 = 1/(g11*g33 - g13*g13); // 1/alpha2
      c2 = g13/g33;
      c3 = g13/g11;
      c4 = 1/g22_y2[iptr2];

      temp1 = -c4*((x_xi-c2*x_zt)*x_etet + (y_xi-c2*y_zt)*y_etet + (z_xi-c2*z_zt)*z_etet);
      temp2 = -c4*(x_et*x_etet + y_et*y_etet + z_et*z_etet);
      temp3 = -c4*((x_zt-c3*x_xi)*x_etet + (y_zt-c3*y_xi)*y_etet + (z_zt-c3*z_xi)*z_etet);
      temp_x = g33*x_xixi + g11*x_ztzt - 2*g13*x_xizt;
      temp_y = g33*y_xixi + g11*y_ztzt - 2*g13*y_xizt;
      temp_z = g33*z_xixi + g11*z_ztzt - 2*g13*z_xizt;

      P[iptr] = temp1 - c1*((x_xi-c2*x_zt)*temp_x + (y_xi-c2*y_zt)*temp_y
                          + (z_xi-c2*z_zt)*temp_z);

      Q[iptr] = temp2 - c1*(x_et*temp_x + y_et*temp_y + z_et*temp_z);

      R[iptr] = temp3 - c1*((x_zt-c3*x_xi)*temp_x + (y_zt-c3*y_xi)*temp_y 
                          + (z_zt-c3*z_xi)*temp_z);
    }
  }
  // bdry z1 
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = 0*siz_iz + j*siz_iy + i;     //(i,j,0)
      iptr1 = 0*siz_iz + j*siz_iy + i+1;   //(i+1,j,0)
      iptr2 = 0*siz_iz + j*siz_iy + i-1;   //(i-1,j,0)
      iptr3 = 0*siz_iz + (j+1)*siz_iy + i; //(i,j+1,0)
      iptr4 = 0*siz_iz + (j-1)*siz_iy + i; //(i,j-1,0)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]); 
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]); 
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]); 
      x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
      y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
      z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
      x_etet = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
      y_etet = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
      z_etet = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

      iptr1 = 0*siz_iz + (j-1)*siz_iy + i-1; //(i-1,j-1,0)
      iptr2 = 0*siz_iz + (j-1)*siz_iy + i+1; //(i+1,j-1,0)
      iptr3 = 0*siz_iz + (j+1)*siz_iy + i+1; //(i+1,j+1,0)
      iptr4 = 0*siz_iz + (j+1)*siz_iy + i-1; //(i-1,j+1,0)

      x_xiet = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
      y_xiet = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
      z_xiet = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
      
      iptr1 = 1*siz_iz + j*siz_iy + i;     //(i,j,1)
      iptr2 = j*nx + i;
      size = nx*ny;
      x_zt = x3d[iptr1] - x3d[iptr];
      y_zt = y3d[iptr1] - y3d[iptr];
      z_zt = z3d[iptr1] - z3d[iptr];
      x_ztzt = p_z1[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
      y_ztzt = p_z1[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
      z_ztzt = p_z1[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

      g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
      g22 = x_et*x_et + y_et*y_et + z_et*z_et;
      g12 = x_xi*x_et + y_xi*y_et + z_xi*z_et;

      c1 = 1/(g11*g22 - g12*g12); // 1/alpha3
      c2 = g12/g22;
      c3 = g12/g11;
      c4 = 1/g33_z1[iptr2];

      temp1 = -c4*((x_xi-c2*x_et)*x_ztzt + (y_xi-c2*y_et)*y_ztzt + (z_xi-c2*z_et)*z_ztzt);
      temp2 = -c4*((x_et-c3*x_xi)*x_ztzt + (y_et-c3*y_xi)*y_ztzt + (z_et-c3*z_xi)*z_ztzt);
      temp3 = -c4*(x_zt*x_ztzt + y_zt*y_ztzt + z_zt*z_ztzt);
      temp_x = g22*x_xixi + g11*x_etet - 2*g12*x_xiet;
      temp_y = g22*y_xixi + g11*y_etet - 2*g12*y_xiet;
      temp_z = g22*z_xixi + g11*z_etet - 2*g12*z_xiet;

      P[iptr] = temp1 - c1*((x_xi-c2*x_et)*temp_x + (y_xi-c2*y_et)*temp_y
                          + (z_xi-c2*z_et)*temp_z);

      Q[iptr] = temp2 - c1*((x_et-c3*x_xi)*temp_x + (y_et-c3*y_xi)*temp_y 
                          + (z_et-c3*z_xi)*temp_z);

      R[iptr] = temp3 - c1*(x_zt*temp_x + y_zt*temp_y + z_zt*temp_z);
    }
  }
  // bdry z2 
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = (nz-1)*siz_iz + j*siz_iy + i;     //(i,j,nz-1)
      iptr1 = (nz-1)*siz_iz + j*siz_iy + i+1;   //(i+1,j,nz-1)
      iptr2 = (nz-1)*siz_iz + j*siz_iy + i-1;   //(i-1,j,nz-1)
      iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i; //(i,j+1,nz-1)
      iptr4 = (nz-1)*siz_iz + (j-1)*siz_iy + i; //(i,j-1,nz-1)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]); 
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]); 
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]); 
      x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
      y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
      z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
      x_etet = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
      y_etet = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
      z_etet = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

      iptr1 = (nz-1)*siz_iz + (j-1)*siz_iy + i-1; //(i-1,j-1,nz-1)
      iptr2 = (nz-1)*siz_iz + (j-1)*siz_iy + i+1; //(i+1,j-1,nz-1)
      iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i+1; //(i+1,j+1,nz-1)
      iptr4 = (nz-1)*siz_iz + (j+1)*siz_iy + i-1; //(i-1,j+1,nz-1)

      x_xiet = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
      y_xiet = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
      z_xiet = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
      
      iptr1 = (nz-2)*siz_iz + j*siz_iy + i;     //(i,j,nz-2)
      iptr2 = j*nx + i;
      size = nx*ny;
      x_zt = x3d[iptr] - x3d[iptr1];
      y_zt = y3d[iptr] - y3d[iptr1];
      z_zt = z3d[iptr] - z3d[iptr1];
      x_ztzt = p_z2[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
      y_ztzt = p_z2[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
      z_ztzt = p_z2[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

      g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
      g22 = x_et*x_et + y_et*y_et + z_et*z_et;
      g12 = x_xi*x_et + y_xi*y_et + z_xi*z_et;

      c1 = 1/(g11*g22 - g12*g12); // 1/alpha3
      c2 = g12/g22;
      c3 = g12/g11;
      c4 = 1/g33_z2[iptr2];

      temp1 = -c4*((x_xi-c2*x_et)*x_ztzt + (y_xi-c2*y_et)*y_ztzt + (z_xi-c2*z_et)*z_ztzt);
      temp2 = -c4*((x_et-c3*x_xi)*x_ztzt + (y_et-c3*y_xi)*y_ztzt + (z_et-c3*z_xi)*z_ztzt);
      temp3 = -c4*(x_zt*x_ztzt + y_zt*y_ztzt + z_zt*z_ztzt);
      temp_x = g22*x_xixi + g11*x_etet - 2*g12*x_xiet;
      temp_y = g22*y_xixi + g11*y_etet - 2*g12*y_xiet;
      temp_z = g22*z_xixi + g11*z_etet - 2*g12*z_xiet;

      P[iptr] = temp1 - c1*((x_xi-c2*x_et)*temp_x + (y_xi-c2*y_et)*temp_y
                          + (z_xi-c2*z_et)*temp_z);

      Q[iptr] = temp2 - c1*((x_et-c3*x_xi)*temp_x + (y_et-c3*y_xi)*temp_y 
                          + (z_et-c3*z_xi)*temp_z);

      R[iptr] = temp3 - c1*(x_zt*temp_x + y_zt*temp_y + z_zt*temp_z);
    }
  }

  interp_inner_source(P,Q,R,nx,ny,nz,coef,first_dire_itype,second_dire_itype);

  return 0;
}

int
ghost_cal(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz, 
          float *p_x, float *p_y, float *p_z,
          float *g11_x, float *g22_y, float *g33_z)
{
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xi0,y_xi0,z_xi0;
  float x_et0,y_et0,z_et0;
  float x_zt0,y_zt0,z_zt0;
  float vn_x,vn_y,vn_z;
  float vn_x0,vn_y0,vn_z0;
  float len_vn,coef;
  float *p_x1,*p_x2;
  float *p_y1,*p_y2;
  float *p_z1,*p_z2;
  float *g11_x1,*g11_x2;
  float *g22_y1,*g22_y2;
  float *g33_z1,*g33_z2;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;
  size_t size;
  p_x1 = p_x;
  p_x2 = p_x + ny*nz*3;
  p_y1 = p_y;
  p_y2 = p_y + nx*nz*3;
  p_z1 = p_z;
  p_z2 = p_z + nx*ny*3;
  g11_x1 = g11_x;
  g11_x2 = g11_x + ny*nz;
  g22_y1 = g22_y;
  g22_y2 = g22_y + nx*nz;
  g33_z1 = g33_z;
  g33_z2 = g33_z + nx*ny;
  // bdry x1 r_et X r_zt
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr1 = k*siz_iz + j*siz_iy + 1;    //(1,j,k)
      iptr2 = k*siz_iz + j*siz_iy + 0;    //(0,j,k)
      x_xi0 = x3d[iptr1] - x3d[iptr2];
      y_xi0 = y3d[iptr1] - y3d[iptr2];
      z_xi0 = z3d[iptr1] - z3d[iptr2];
      iptr1 = k*siz_iz + (j+1)*siz_iy + 0;   //(0,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + 0;   //(0,j-1,k)
      x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
      iptr3 = (k+1)*siz_iz + j*siz_iy + 0;   //(0,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + 0;   //(0,j,k-1)
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
      // orth vector. cross product
      vn_x = y_et*z_zt - z_et*y_zt;
      vn_y = z_et*x_zt - x_et*z_zt;
      vn_z = x_et*y_zt - y_et*x_zt;
      len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      // norm
      vn_x0 = vn_x/len_vn;
      vn_y0 = vn_y/len_vn;
      vn_z0 = vn_z/len_vn;
      // projection from r_xi0 to vn. dot product
      coef = x_xi0*vn_x0 + y_xi0*vn_y0 + z_xi0*vn_z0;

      x_xi = coef*vn_x0;
      y_xi = coef*vn_y0;
      z_xi = coef*vn_z0;
       
      iptr = k*ny + j;
      iptr1 = k*siz_iz + j*siz_iy + 0;
      size = ny*nz;
      p_x1[iptr+0*size] = x3d[iptr1] - x_xi;
      p_x1[iptr+1*size] = y3d[iptr1] - y_xi;
      p_x1[iptr+2*size] = z3d[iptr1] - z_xi;
      g11_x1[iptr] = pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2);
    }
  }
  // bdry x2 r_et X r_zt
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr1 = k*siz_iz + j*siz_iy + nx-1;    //(nx-1,j,k)
      iptr2 = k*siz_iz + j*siz_iy + nx-2;    //(nx-2,j,k)
      x_xi0 = x3d[iptr1] - x3d[iptr2];
      y_xi0 = y3d[iptr1] - y3d[iptr2];
      z_xi0 = z3d[iptr1] - z3d[iptr2];
      iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1;   //(nx-1,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1;   //(nx-1,j-1,k)
      x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
      iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1;   //(nx-1,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + nx-1;   //(nx-1,j,k-1)
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
      // orth vector. cross product
      vn_x = y_et*z_zt - z_et*y_zt;
      vn_y = z_et*x_zt - x_et*z_zt;
      vn_z = x_et*y_zt - y_et*x_zt;
      len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      // norm
      vn_x0 = vn_x/len_vn;
      vn_y0 = vn_y/len_vn;
      vn_z0 = vn_z/len_vn;
      // projection from r_xi0 to vn. dot product
      coef = x_xi0*vn_x0 + y_xi0*vn_y0 + z_xi0*vn_z0;

      x_xi = coef*vn_x0;
      y_xi = coef*vn_y0;
      z_xi = coef*vn_z0;
       
      iptr = k*ny + j;
      iptr1 = k*siz_iz + j*siz_iy + nx-1;
      size = ny*nz;
      p_x2[iptr+0*size] = x3d[iptr1] + x_xi;
      p_x2[iptr+1*size] = y3d[iptr1] + y_xi;
      p_x2[iptr+2*size] = z3d[iptr1] + z_xi;
      g11_x2[iptr] = pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2);
    }
  }
  // bdry y1 r_zt X r_xi
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = k*siz_iz + 1*siz_iy + i;    //(i,1,k)
      iptr2 = k*siz_iz + 0*siz_iy + i;    //(i,0,k)
      x_et0 = x3d[iptr1] - x3d[iptr2];
      y_et0 = y3d[iptr1] - y3d[iptr2];
      z_et0 = z3d[iptr1] - z3d[iptr2];
      iptr1 = k*siz_iz + 0*siz_iy + (i+1);   //(i+1,0,k)
      iptr2 = k*siz_iz + 0*siz_iy + (i-1);   //(i-1,0,k)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
      iptr3 = (k+1)*siz_iz + 0*siz_iy + i;   //(i,0,k+1)
      iptr4 = (k-1)*siz_iz + 0*siz_iy + i;   //(i,0,k-1)
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
      // orth vector. cross product
      vn_x = y_zt*z_xi - z_zt*y_xi;
      vn_y = z_zt*x_xi - x_zt*z_xi;
      vn_z = x_zt*y_xi - y_zt*x_xi;
      len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      // norm
      vn_x0 = vn_x/len_vn;
      vn_y0 = vn_y/len_vn;
      vn_z0 = vn_z/len_vn;
      // projection from r_et0 to vn. dot product
      coef = x_et0*vn_x0 + y_et0*vn_y0 + z_et0*vn_z0;

      x_et = coef*vn_x0;
      y_et = coef*vn_y0;
      z_et = coef*vn_z0;
       
      iptr = k*nx + i;
      iptr1 = k*siz_iz + 0*siz_iy + i;
      size =  nx*nz;
      p_y1[iptr+0*size] = x3d[iptr1] - x_et;
      p_y1[iptr+1*size] = y3d[iptr1] - y_et;
      p_y1[iptr+2*size] = z3d[iptr1] - z_et;
      g22_y1[iptr] = pow(x_et,2) + pow(y_et,2) + pow(z_et,2);
    }
  }
  // bdry y2 r_zt X r_xi
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = k*siz_iz + (ny-1)*siz_iy + i;    //(i,ny-1,k)
      iptr2 = k*siz_iz + (ny-2)*siz_iy + i;    //(i,ny-2,k)
      x_et0 = x3d[iptr1] - x3d[iptr2];
      y_et0 = y3d[iptr1] - y3d[iptr2];
      z_et0 = z3d[iptr1] - z3d[iptr2];
      iptr1 = k*siz_iz + (ny-1)*siz_iy + (i+1);   //(i+1,ny-1,k)
      iptr2 = k*siz_iz + (ny-1)*siz_iy + (i-1);   //(i-1,ny-1,k)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
      iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k+1)
      iptr4 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k-1)
      x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
      // orth vector. cross product
      vn_x = y_zt*z_xi - z_zt*y_xi;
      vn_y = z_zt*x_xi - x_zt*z_xi;
      vn_z = x_zt*y_xi - y_zt*x_xi;
      len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      // norm
      vn_x0 = vn_x/len_vn;
      vn_y0 = vn_y/len_vn;
      vn_z0 = vn_z/len_vn;
      // projection from r_et0 to vn. dot product
      coef = x_et0*vn_x0 + y_et0*vn_y0 + z_et0*vn_z0;

      x_et = coef*vn_x0;
      y_et = coef*vn_y0;
      z_et = coef*vn_z0;
       
      iptr = k*nx + i;
      iptr1 = k*siz_iz + (ny-1)*siz_iy + i;
      size = nx*nz;
      p_y2[iptr+0*size] = x3d[iptr1] + x_et;
      p_y2[iptr+1*size] = y3d[iptr1] + y_et;
      p_y2[iptr+2*size] = z3d[iptr1] + z_et;
      g22_y2[iptr] = pow(x_et,2) + pow(y_et,2) + pow(z_et,2);
    }
  }
  // bdry z1 r_xi X r_et
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = 1*siz_iz + j*siz_iy + i;    //(i,j,1)
      iptr2 = 0*siz_iz + j*siz_iy + i;    //(i,j,0)
      x_zt0 = x3d[iptr1] - x3d[iptr2];
      y_zt0 = y3d[iptr1] - y3d[iptr2];
      z_zt0 = z3d[iptr1] - z3d[iptr2];
      iptr1 = 0*siz_iz + j*siz_iy + (i+1);   //(i+1,j,0)
      iptr2 = 0*siz_iz + j*siz_iy + (i-1);   //(i-1,j,0)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
      iptr3 = 0*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,0)
      iptr4 = 0*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,0)
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);
      // orth vector. cross product
      vn_x = y_xi*z_et - z_xi*y_et;
      vn_y = z_xi*x_et - x_xi*z_et;
      vn_z = x_xi*y_et - y_xi*x_et;
      len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      // norm
      vn_x0 = vn_x/len_vn;
      vn_y0 = vn_y/len_vn;
      vn_z0 = vn_z/len_vn;
      // projection from r_et0 to vn. dot product
      coef = x_zt0*vn_x0 + y_zt0*vn_y0 + z_zt0*vn_z0;

      x_zt = coef*vn_x0;
      y_zt = coef*vn_y0;
      z_zt = coef*vn_z0;
       
      iptr = j*nx + i;
      iptr1 = 0*siz_iz + j*siz_iy + i;
      size = nx*ny;
      p_z1[iptr+0*size] = x3d[iptr1] - x_zt;
      p_z1[iptr+1*size] = y3d[iptr1] - y_zt;
      p_z1[iptr+2*size] = z3d[iptr1] - z_zt;
      g33_z1[iptr] = pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2);
    }
  }
  // bdry z2 r_xi X r_et
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (nz-1)*siz_iz + j*siz_iy + i;    //(i,j,nz-1)
      iptr2 = (nz-2)*siz_iz + j*siz_iy + i;    //(i,j,nz-2)
      x_zt0 = x3d[iptr1] - x3d[iptr2];
      y_zt0 = y3d[iptr1] - y3d[iptr2];
      z_zt0 = z3d[iptr1] - z3d[iptr2];
      iptr1 = (nz-1)*siz_iz + j*siz_iy + (i+1);   //(i+1,j,nz-1)
      iptr2 = (nz-1)*siz_iz + j*siz_iy + (i-1);   //(i-1,j,nz-1)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
      iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,nz-1)
      iptr4 = (nz-1)*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,nz-1)
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);
      // orth vector. cross product
      vn_x = y_xi*z_et - z_xi*y_et;
      vn_y = z_xi*x_et - x_xi*z_et;
      vn_z = x_xi*y_et - y_xi*x_et;
      len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      // norm
      vn_x0 = vn_x/len_vn;
      vn_y0 = vn_y/len_vn;
      vn_z0 = vn_z/len_vn;
      // projection from r_et0 to vn. dot product
      coef = x_zt0*vn_x0 + y_zt0*vn_y0 + z_zt0*vn_z0;

      x_zt = coef*vn_x0;
      y_zt = coef*vn_y0;
      z_zt = coef*vn_z0;
       
      iptr = j*nx + i;
      iptr1 = (nz-1)*siz_iz + j*siz_iy + i;
      size = nx*ny;
      p_z2[iptr+0*size] = x3d[iptr1] + x_zt;
      p_z2[iptr+1*size] = y3d[iptr1] + y_zt;
      p_z2[iptr+2*size] = z3d[iptr1] + z_zt;
      g33_z2[iptr] = pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2);
    }
  }

  return 0;
}

int
higen_gene(gd_t *gdcurv, par_t *par)
{
  float err = par->iter_err;
  int max_iter = par->max_iter;
  float coef = par->coef;
  int first_dire_itype = par->first_dire_itype;
  int second_dire_itype = par->second_dire_itype;

  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  float *P = NULL; //source term P
  float *Q = NULL; //source term Q
  float *R = NULL; //source term Q
  P = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source P");
  Q = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source Q");
  R = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source R");

  float *x3d_tmp =NULL;
  float *y3d_tmp =NULL;
  float *z3d_tmp =NULL;
  x3d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "x3d_tmp");
  y3d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "y3d_tmp");
  z3d_tmp = (float *)mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "z3d_tmp");

  float dx1,dx2,dy1,dy2,dz1,dz2;
  dx1 = par->distance[0];
  dx2 = par->distance[1];
  dy1 = par->distance[2];
  dy2 = par->distance[3];
  dz1 = par->distance[4];
  dz2 = par->distance[5];
  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  size_t iptr, iptr1, iptr2, iptr3;
  float resi, resj, resk;
  float max_resi, max_resj, max_resk;
  
  // copy coord
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + j*siz_iy + i;
        x3d_tmp[iptr] = x3d[iptr];
        y3d_tmp[iptr] = y3d[iptr];
        z3d_tmp[iptr] = z3d[iptr];
      }
    }
  }

  int flag_true = 1;
  while(flag_true)
  {
    // update solver
    update_SOR(x3d,y3d,z3d,x3d_tmp,y3d_tmp,z3d_tmp,nx,ny,nz,P,Q,R,w);
    Niter += 1;

    max_resi = 0.0;
    max_resj = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + i+1;
          iptr2 = k*siz_iz + (j+1)*siz_iy + i;
          iptr3 = (k+1)*siz_iz + j*siz_iy + i;
          resi = fabs((x3d_tmp[iptr] - x3d[iptr])/(x3d[iptr1]-x3d[iptr]));
          resj = fabs((y3d_tmp[iptr] - y3d[iptr])/(y3d[iptr2]-y3d[iptr]));
          resk = fabs((z3d_tmp[iptr] - z3d[iptr])/(z3d[iptr3]-z3d[iptr]));
          max_resi = fmax(max_resi,resi);
          max_resj = fmax(max_resj,resj);
          max_resk = fmax(max_resk,resk);
        }
      }
    }

    // copy coord
    for(int k=0; k<nz; k++) {
      for(int j=0; j<ny; j++) {
        for(int i=0; i<nx; i++)
        {
          iptr = k*siz_iz + j*siz_iy + i;
          x3d[iptr] = x3d_tmp[iptr];
          y3d[iptr] = y3d_tmp[iptr];
          z3d[iptr] = z3d_tmp[iptr];
        }
      }
    }
    
    set_src_higen(x3d,y3d,z3d,nx,ny,nz,P,Q,R,dx1,dx2,dy1,dy2,dz1,dz2,
                  coef,first_dire_itype,second_dire_itype);

    fprintf(stdout,"number of iter is %d\n", Niter);
    fprintf(stdout,"max_resi is %f, max_resj is %f, max_resk is %f\n", max_resi, max_resj, max_resk);
    fflush(stdout);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resj < err && max_resk < err) {
      flag_true = 0;
    }

  }

  free(x3d_tmp);
  free(y3d_tmp);
  free(z3d_tmp);
  free(P);
  free(Q);
  free(R);

  return 0;
}

int
set_src_higen(float *x3d, float *y3d, float *z3d, int nx, int ny, int nz,
              float *P, float *Q, float *R, float dx1, float dx2, 
              float dy1, float dy2, float dz1, float dz2,
              float coef, int first_dire_itype, int second_dire_itype)
{
  float theta0 = PI/2;
  float a = 0.1;
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float dot, len_xi, len_et, len_zt, dif_dis;
  float cos_theta, theta, dif_theta;
  float x_xi, y_xi, z_xi;
  float x_et, y_et, z_et;
  float x_zt, y_zt, z_zt;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;

  // NOTE: P Q R update need change sign
  // when xi = 1, et = 1, zt = 1. 
  // for example:
  // xi = 0. Q R -, P + 
  // xi = 1. Q R +, P - 
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr  = k*siz_iz + j*siz_iy + 0;         //(0,j,k)
      iptr1 = k*siz_iz + (j+1)*siz_iy + 0;     //(0,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + 0;     //(0,j-1,k)
      iptr3 = (k+1)*siz_iz + j*siz_iy + 0;     //(0,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + 0;     //(0,j,k-1)

      x_et = 0.5*(x3d[iptr1]-x3d[iptr2]);
      y_et = 0.5*(y3d[iptr1]-y3d[iptr2]);
      z_et = 0.5*(z3d[iptr1]-z3d[iptr2]);
      x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

      iptr1 = k*siz_iz + j*siz_iy + 1;     //(1,j,k)

      x_xi = x3d[iptr1] - x3d[iptr];
      y_xi = y3d[iptr1] - y3d[iptr];
      z_xi = z3d[iptr1] - z3d[iptr];
      
      // cos(theta) = a.b/(|a|*|b|)
      len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
      len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
      len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
      // r_xi . r_et
      dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
      cos_theta = dot/(len_xi*len_et);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      Q[iptr] = Q[iptr] - a*tanh(dif_theta);

      // r_xi . r_zt
      dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      R[iptr] = R[iptr] - a*tanh(dif_theta);

      dif_dis = (dx1-len_xi)/dx1;
      P[iptr] = P[iptr] + a*tanh(dif_dis);
    }
  }
  // bdry x2 xi=1
  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++)
    {
      iptr  = k*siz_iz + j*siz_iy + nx-1;         //(nx-1,j,k)
      iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1;     //(nx-1,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1;     //(nx-1,j-1,k)
      iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1;     //(nx-1,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + nx-1;     //(nx-1,j,k-1)

      x_et = 0.5*(x3d[iptr1]-x3d[iptr2]);
      y_et = 0.5*(y3d[iptr1]-y3d[iptr2]);
      z_et = 0.5*(z3d[iptr1]-z3d[iptr2]);
      x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

      iptr1 = k*siz_iz + j*siz_iy + nx-2;     //(nx-2,j,k)

      x_xi = x3d[iptr] - x3d[iptr1];
      y_xi = y3d[iptr] - y3d[iptr1];
      z_xi = z3d[iptr] - z3d[iptr1];
      
      // cos(theta) = a.b/(|a|*|b|)
      len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
      len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
      len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
      // r_xi . r_et
      dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
      cos_theta = dot/(len_xi*len_et);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      Q[iptr] = Q[iptr] + a*tanh(dif_theta);

      // r_xi . r_zt
      dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      R[iptr] = R[iptr] + a*tanh(dif_theta);

      dif_dis = (dx2-len_xi)/dx2;
      P[iptr] = P[iptr] - a*tanh(dif_dis);
    }
  }
  // bdry y1 et=0
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = k*siz_iz + 0*siz_iy + i;       //(i,0,k)
      iptr1 = k*siz_iz + 0*siz_iy + i+1;     //(i+1,0,k)
      iptr2 = k*siz_iz + 0*siz_iy + i-1;     //(i-1,0,k)
      iptr3 = (k+1)*siz_iz + 0*siz_iy + i;   //(i,0,k+1)
      iptr4 = (k-1)*siz_iz + 0*siz_iy + i;   //(i,0,k-1)

      x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
      x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

      iptr1 = k*siz_iz + 1*siz_iy + i;     //(i,1,k)

      x_et = x3d[iptr1] - x3d[iptr];
      y_et = y3d[iptr1] - y3d[iptr];
      z_et = z3d[iptr1] - z3d[iptr];
      
      // cos(theta) = a.b/(|a|*|b|)
      len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
      len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
      len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
      // r_xi . r_et
      dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
      cos_theta = dot/(len_xi*len_et);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      P[iptr] = P[iptr] - a*tanh(dif_theta);

      // r_et . r_zt
      dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
      cos_theta = dot/(len_et*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      R[iptr] = R[iptr] - a*tanh(dif_theta);

      dif_dis = (dy1-len_et)/dy1;
      Q[iptr] = Q[iptr] + a*tanh(dif_dis);
    }
  }
  // bdry y2 et=0
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = k*siz_iz + (ny-1)*siz_iy + i;       //(i,ny-1,k)
      iptr1 = k*siz_iz + (ny-1)*siz_iy + i+1;     //(i+1,ny-1,k)
      iptr2 = k*siz_iz + (ny-1)*siz_iy + i-1;     //(i-1,ny-1,k)
      iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k+1)
      iptr4 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k-1)

      x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
      x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
      y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
      z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

      iptr1 = k*siz_iz + (ny-2)*siz_iy + i;     //(i,ny-2,k)

      x_et = x3d[iptr] - x3d[iptr1];
      y_et = y3d[iptr] - y3d[iptr1];
      z_et = z3d[iptr] - z3d[iptr1];
      
      // cos(theta) = a.b/(|a|*|b|)
      len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
      len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
      len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
      // r_xi . r_et
      dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
      cos_theta = dot/(len_xi*len_et);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      P[iptr] = P[iptr] + a*tanh(dif_theta);

      // r_et . r_zt
      dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
      cos_theta = dot/(len_et*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      R[iptr] = R[iptr] + a*tanh(dif_theta);

      dif_dis = (dy2-len_et)/dy2;
      Q[iptr] = Q[iptr] - a*tanh(dif_dis);
    }
  }
  // bdry z1 zt=0
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = 0*siz_iz + j*siz_iy + i;       //(i,j,0)
      iptr1 = 0*siz_iz + j*siz_iy + i+1;     //(i+1,j,0)
      iptr2 = 0*siz_iz + j*siz_iy + i-1;     //(i-1,j,0)
      iptr3 = 0*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,0)
      iptr4 = 0*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,0)

      x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
      x_et = 0.5*(x3d[iptr3]-x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3]-y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3]-z3d[iptr4]);

      iptr1 = 1*siz_iz + j*siz_iy + i;     //(i,j,1)

      x_zt = x3d[iptr1] - x3d[iptr];
      y_zt = y3d[iptr1] - y3d[iptr];
      z_zt = z3d[iptr1] - z3d[iptr];
      
      // cos(theta) = a.b/(|a|*|b|)
      len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
      len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
      len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
      // r_xi . r_zt
      dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      P[iptr] = P[iptr] - a*tanh(dif_theta);

      // r_et . r_zt
      dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
      cos_theta = dot/(len_et*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      Q[iptr] = Q[iptr] - a*tanh(dif_theta);

      dif_dis = (dz1-len_zt)/dz1;
      R[iptr] = R[iptr] + a*tanh(dif_dis);
    }
  }
  // bdry z2 zt=1
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = (nz-1)*siz_iz + j*siz_iy + i;       //(i,j,nz-1)
      iptr1 = (nz-1)*siz_iz + j*siz_iy + i+1;     //(i+1,j,nz-1)
      iptr2 = (nz-1)*siz_iz + j*siz_iy + i-1;     //(i-1,j,nz-1)
      iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,nz-1)
      iptr4 = (nz-1)*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,nz-1)

      x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
      x_et = 0.5*(x3d[iptr3]-x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3]-y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3]-z3d[iptr4]);

      iptr1 = (nz-2)*siz_iz + j*siz_iy + i;     //(i,j,nz-2)

      x_zt = x3d[iptr] - x3d[iptr1];
      y_zt = y3d[iptr] - y3d[iptr1];
      z_zt = z3d[iptr] - z3d[iptr1];
      
      // cos(theta) = a.b/(|a|*|b|)
      len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
      len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
      len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
      // r_xi . r_zt
      dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      P[iptr] = P[iptr] + a*tanh(dif_theta);

      // r_et . r_zt
      dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
      cos_theta = dot/(len_et*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      Q[iptr] = Q[iptr] + a*tanh(dif_theta);

      dif_dis = (dz2-len_zt)/dz2;
      R[iptr] = R[iptr] - a*tanh(dif_dis);
    }
  }

  interp_inner_source(P,Q,R,nx,ny,nz,coef,first_dire_itype,second_dire_itype);

  return 0;
}

int
interp_inner_source(float *P, float *Q, float *R,
                    int nx, int ny, int nz, float coef,
                    int first_dire_itype,int second_dire_itype)
{
  size_t siz_iz = nx*ny;
  size_t siz_iy = nx;
  float xi,et,zt,c0,c1,r0,r1;
  size_t iptr,iptr1,iptr2;
  if(first_dire_itype == Z_DIRE && second_dire_itype == Y_DIRE)
  {
    // z direciton
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          zt = (1.0*k)/(nz-1);
          c0 = 1-zt;
          c1 = zt;

          r0 = exp(coef*zt);
          r1 = exp(coef*(1-zt));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = 0*siz_iz + j*siz_iy + i;  
          iptr2 = (nz-1)*siz_iz + j*siz_iy + i;  

          P[iptr] = r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = c0*R[iptr1] + c1*R[iptr2];
        }
      }
    }
    // y direction
    for(int k=5; k<nz-5; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          et = (1.0*j)/(ny-1);
          c0 = 1-et;
          c1 = et;

          r0 = exp(coef*et);
          r1 = exp(coef*(1-et));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + 0*siz_iy + i;  
          iptr2 = k*siz_iz + (ny-1)*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + c0*Q[iptr1] + c1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // x direction
    for(int k=5; k<nz-5; k++) {
      for(int j=5; j<ny-5; j++) {
        for(int i=1; i<nx-1; i++)
        {
          xi = (1.0*i)/(nx-1);
          c0 = 1-xi;
          c1 = xi;

          r0 = exp(coef*xi);
          r1 = exp(coef*(1-xi));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + 0;  
          iptr2 = k*siz_iz + j*siz_iy + nx-1;  

          P[iptr] = P[iptr] + c0*P[iptr1] + c1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
  }
  if(first_dire_itype == Z_DIRE && second_dire_itype == X_DIRE)
  {
    // z direciton
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          zt = (1.0*k)/(nz-1);
          c0 = 1-zt;
          c1 = zt;

          r0 = exp(coef*zt);
          r1 = exp(coef*(1-zt));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = 0*siz_iz + j*siz_iy + i;  
          iptr2 = (nz-1)*siz_iz + j*siz_iy + i;  

          P[iptr] = r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = c0*R[iptr1] + c1*R[iptr2];
        }
      }
    }
    // x direction
    for(int k=5; k<nz-5; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          xi = (1.0*i)/(nx-1);
          c0 = 1-xi;
          c1 = xi;

          r0 = exp(coef*xi);
          r1 = exp(coef*(1-xi));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + 0;  
          iptr2 = k*siz_iz + j*siz_iy + nx-1;  

          P[iptr] = P[iptr] + c0*P[iptr1] + c1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // y direction
    for(int k=5; k<nz-5; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=5; i<nx-5; i++)
        {
          et = (1.0*j)/(ny-1);
          c0 = 1-et;
          c1 = et;

          r0 = exp(coef*et);
          r1 = exp(coef*(1-et));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + 0*siz_iy + i;  
          iptr2 = k*siz_iz + (ny-1)*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + c0*Q[iptr1] + c1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
  }
  if(first_dire_itype == Y_DIRE && second_dire_itype == Z_DIRE)
  {
    // y direction
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          et = (1.0*j)/(ny-1);
          c0 = 1-et;
          c1 = et;

          r0 = exp(coef*et);
          r1 = exp(coef*(1-et));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + 0*siz_iy + i;  
          iptr2 = k*siz_iz + (ny-1)*siz_iy + i;  

          P[iptr] = r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = c0*Q[iptr1] + c1*Q[iptr2];
          R[iptr] = r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // z direciton
    for(int k=1; k<nz-1; k++) {
      for(int j=5; j<ny-5; j++) {
        for(int i=1; i<nx-1; i++)
        {
          zt = (1.0*k)/(nz-1);
          c0 = 1-zt;
          c1 = zt;

          r0 = exp(coef*zt);
          r1 = exp(coef*(1-zt));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = 0*siz_iz + j*siz_iy + i;  
          iptr2 = (nz-1)*siz_iz + j*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + c0*R[iptr1] + c1*R[iptr2];
        }
      }
    }
    // x direction
    for(int k=5; k<nz-5; k++) {
      for(int j=5; j<ny-5; j++) {
        for(int i=1; i<nx-1; i++)
        {
          xi = (1.0*i)/(nx-1);
          c0 = 1-xi;
          c1 = xi;

          r0 = exp(coef*xi);
          r1 = exp(coef*(1-xi));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + 0;  
          iptr2 = k*siz_iz + j*siz_iy + nx-1;  

          P[iptr] = P[iptr] + c0*P[iptr1] + c1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
  }
  if(first_dire_itype == Y_DIRE && second_dire_itype == X_DIRE)
  {
    // y direction
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          et = (1.0*j)/(ny-1);
          c0 = 1-et;
          c1 = et;

          r0 = exp(coef*et);
          r1 = exp(coef*(1-et));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + 0*siz_iy + i;  
          iptr2 = k*siz_iz + (ny-1)*siz_iy + i;  

          P[iptr] = r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = c0*Q[iptr1] + c1*Q[iptr2];
          R[iptr] = r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // x direction
    for(int k=1; k<nz-1; k++) {
      for(int j=5; j<ny-5; j++) {
        for(int i=1; i<nx-1; i++)
        {
          xi = (1.0*i)/(nx-1);
          c0 = 1-xi;
          c1 = xi;

          r0 = exp(coef*xi);
          r1 = exp(coef*(1-xi));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + 0;  
          iptr2 = k*siz_iz + j*siz_iy + nx-1;  

          P[iptr] = P[iptr] + c0*P[iptr1] + c1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // z direciton
    for(int k=1; k<nz-1; k++) {
      for(int j=5; j<ny-5; j++) {
        for(int i=5; i<nx-5; i++)
        {
          zt = (1.0*k)/(nz-1);
          c0 = 1-zt;
          c1 = zt;

          r0 = exp(coef*zt);
          r1 = exp(coef*(1-zt));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = 0*siz_iz + j*siz_iy + i;  
          iptr2 = (nz-1)*siz_iz + j*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + c0*R[iptr1] + c1*R[iptr2];
        }
      }
    }
  }
  if(first_dire_itype == X_DIRE && second_dire_itype == Y_DIRE)
  {
    // x direction
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          xi = (1.0*i)/(nx-1);
          c0 = 1-xi;
          c1 = xi;

          r0 = exp(coef*xi);
          r1 = exp(coef*(1-xi));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + 0;  
          iptr2 = k*siz_iz + j*siz_iy + nx-1;  

          P[iptr] = c0*P[iptr1] + c1*P[iptr2];
          Q[iptr] = r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // y direction
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=5; i<nx-5; i++)
        {
          et = (1.0*j)/(ny-1);
          c0 = 1-et;
          c1 = et;

          r0 = exp(coef*et);
          r1 = exp(coef*(1-et));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + 0*siz_iy + i;  
          iptr2 = k*siz_iz + (ny-1)*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + c0*Q[iptr1] + c1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // z direciton
    for(int k=1; k<nz-1; k++) {
      for(int j=5; j<ny-5; j++) {
        for(int i=5; i<nx-5; i++)
        {
          zt = (1.0*k)/(nz-1);
          c0 = 1-zt;
          c1 = zt;

          r0 = exp(coef*zt);
          r1 = exp(coef*(1-zt));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = 0*siz_iz + j*siz_iy + i;  
          iptr2 = (nz-1)*siz_iz + j*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + c0*R[iptr1] + c1*R[iptr2];
        }
      }
    }
  }
  if(first_dire_itype == X_DIRE && second_dire_itype == Z_DIRE)
  {
    // x direction
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          xi = (1.0*i)/(nx-1);
          c0 = 1-xi;
          c1 = xi;

          r0 = exp(coef*xi);
          r1 = exp(coef*(1-xi));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + 0;  
          iptr2 = k*siz_iz + j*siz_iy + nx-1;  

          P[iptr] = c0*P[iptr1] + c1*P[iptr2];
          Q[iptr] = r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
    // z direciton
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=5; i<nx-5; i++)
        {
          zt = (1.0*k)/(nz-1);
          c0 = 1-zt;
          c1 = zt;

          r0 = exp(coef*zt);
          r1 = exp(coef*(1-zt));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = 0*siz_iz + j*siz_iy + i;  
          iptr2 = (nz-1)*siz_iz + j*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
          R[iptr] = R[iptr] + c0*R[iptr1] + c1*R[iptr2];
        }
      }
    }
    // y direction
    for(int k=5; k<nz-5; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=5; i<nx-5; i++)
        {
          et = (1.0*j)/(ny-1);
          c0 = 1-et;
          c1 = et;

          r0 = exp(coef*et);
          r1 = exp(coef*(1-et));

          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + 0*siz_iy + i;  
          iptr2 = k*siz_iz + (ny-1)*siz_iy + i;  

          P[iptr] = P[iptr] + r0*P[iptr1] + r1*P[iptr2];
          Q[iptr] = Q[iptr] + c0*Q[iptr1] + c1*Q[iptr2];
          R[iptr] = R[iptr] + r0*R[iptr1] + r1*R[iptr2];
        }
      }
    }
  }

  return 0;
}
