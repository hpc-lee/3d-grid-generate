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
int linear_tfi(gd_t *gdcurv, bdry_t *bdry, mympi_t *mympi)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int gni1 = gdcurv->gni1;
  int gnj1 = gdcurv->gnj1;
  int gnk1 = gdcurv->gnk1;
  int gni2 = gdcurv->gni2;
  int gnj2 = gdcurv->gnj2;
  int gnk2 = gdcurv->gnk2;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  int total_nx = bdry->number_of_grid_points_x;
  int total_ny = bdry->number_of_grid_points_y;
  int total_nz = bdry->number_of_grid_points_z;
  int siz_bx = total_ny*total_nz;
  int siz_by = total_nx*total_nz;
  int siz_bz = total_nx*total_ny;
  float *x1 = bdry->x1;
  float *x2 = bdry->x2;
  float *y1 = bdry->y1;
  float *y2 = bdry->y2;
  float *z1 = bdry->z1;
  float *z2 = bdry->z2;

  int gni, gnj, gnk;
  float xi,et,zt;
  float a0,a1,b0,b1,c0,c1;
  size_t iptr, iptr1, iptr2, iptr3, iptr4;
  size_t iptr5, iptr6, iptr7, iptr8, iptr9;
  size_t iptr10, iptr11, iptr12;
  float U_x,V_x,W_x,UW_x,UV_x,VW_x,UVW_x;
  float U_y,V_y,W_y,UW_y,UV_y,VW_y,UVW_y;
  float U_z,V_z,W_z,UW_z,UV_z,VW_z,UVW_z;
  
  for (int k=nk1; k<=nk2; k++) {
    for (int j=nj1; j<=nj2; j++) {
      for (int i=ni1; i<=ni2; i++)
      {
        gni = gni1 + i; 
        gnj = gnj1 + j; 
        gnk = gnk1 + k; 
        xi = (1.0*gni)/(total_nx-1); // 1.0* equal int to float
        et = (1.0*gnj)/(total_ny-1); 
        zt = (1.0*gnk)/(total_nz-1);
        a0 = 1-xi;
        a1 = xi;
        b0 = 1-et;
        b1 = et;
        c0 = 1-zt;
        c1 = zt;
  
  
        // 6 face points
        iptr1 = gnk*total_ny + gnj;
        iptr2 = gnk*total_nx + gni;
        iptr3 = gnj*total_nx + gni;
        U_x = a0*x1[iptr1] + a1*x2[iptr1];
        V_x = b0*y1[iptr2] + b1*y2[iptr2];
        W_x = c0*z1[iptr3] + c1*z2[iptr3];
  
        U_y = a0*x1[iptr1+siz_bx] + a1*x2[iptr1+siz_bx];
        V_y = b0*y1[iptr2+siz_by] + b1*y2[iptr2+siz_by];
        W_y = c0*z1[iptr3+siz_bz] + c1*z2[iptr3+siz_bz];
  
        U_z = a0*x1[iptr1+2*siz_bx] + a1*x2[iptr1+2*siz_bx];
        V_z = b0*y1[iptr2+2*siz_by] + b1*y2[iptr2+2*siz_by];
        W_z = c0*z1[iptr3+2*siz_bz] + c1*z2[iptr3+2*siz_bz];

        // 12 edge points
  
        iptr1 = gnk*total_ny + 0;  //x1
        iptr2 = gnk*total_ny + total_ny-1; //x1
        iptr3 = gnk*total_ny + 0; //x2
        iptr4 = gnk*total_ny + total_ny-1; //x2
        UV_x = a0*b0*x1[iptr1] + a0*b1*x1[iptr2]
             + a1*b0*x2[iptr3] + a1*b1*x2[iptr4];
        UV_y = a0*b0*x1[iptr1+siz_bx] + a0*b1*x1[iptr2+siz_bx]
             + a1*b0*x2[iptr3+siz_bx] + a1*b1*x2[iptr4+siz_bx];
        UV_z = a0*b0*x1[iptr1+2*siz_bx] + a0*b1*x1[iptr2+2*siz_bx]
             + a1*b0*x2[iptr3+2*siz_bx] + a1*b1*x2[iptr4+2*siz_bx];

        iptr5 = gnj*total_nx + 0;  // z1
        iptr6 = gnj*total_nx + total_nx-1; // z1
        iptr7 = gnj*total_nx + 0;  // z2
        iptr8 = gnj*total_nx + total_nx-1;  // z2
        UW_x = a0*c0*z1[iptr5] + a1*c0*z1[iptr6]
             + a0*c1*z2[iptr7] + a1*c1*z2[iptr8];
        UW_y = a0*c0*z1[iptr5+siz_bz] + a1*c0*z1[iptr6+siz_bz]
             + a0*c1*z2[iptr7+siz_bz] + a1*c1*z2[iptr8+siz_bz];
        UW_z = a0*c0*z1[iptr5+2*siz_bz] + a1*c0*z1[iptr6+2*siz_bz]
             + a0*c1*z2[iptr7+2*siz_bz] + a1*c1*z2[iptr8+2*siz_bz];

        iptr9  = gni; // y1
        iptr10 = (total_nz-1)*total_nx + gni; // y1
        iptr11 = gni; // y2
        iptr12 = (total_nz-1)*total_nx + gni; // y2 
        VW_x = b0*c0*y1[iptr9]  + b0*c1*y1[iptr10]
             + b1*c0*y2[iptr11] + b1*c1*y2[iptr12];
  
        VW_y = b0*c0*y1[iptr9 +siz_by] + b0*c1*y1[iptr10+siz_by]
             + b1*c0*y2[iptr11+siz_by] + b1*c1*y2[iptr12+siz_by];
  
        VW_z = b0*c0*y1[iptr9 +2*siz_by] + b0*c1*y1[iptr10+2*siz_by]
             + b1*c0*y2[iptr11+2*siz_by] + b1*c1*y2[iptr12+2*siz_by];
  
        // 8 corner point
        iptr1 = 0; // x1 (0,0,0)
        iptr2 = 0; // x2 (0,0,nx-1)
        iptr3 = total_ny-1; // x2 (0,ny-1,nx-1)
        iptr4 = total_ny-1; // x1 (0,ny-1,0)
        iptr5 = (total_nz-1)*total_ny; //x1 (nz-1,0,0)
        iptr6 = (total_nz-1)*total_ny; //x2 (nz-1,0,nx-1)
        iptr7 = (total_nz-1)*total_ny + total_ny-1; //x2 (nz-1,ny-1,nx-1)
        iptr8 = (total_nz-1)*total_ny + total_ny-1; //x1 (nz-1,ny-1,0)
        UVW_x = a0*b0*c0*x1[iptr1] + a1*b0*c0*x2[iptr2]
               +a1*b1*c0*x2[iptr3] + a0*b1*c0*x1[iptr4]
               +a0*b0*c1*x1[iptr5] + a1*b0*c1*x2[iptr6]
               +a1*b1*c1*x2[iptr7] + a0*b1*c1*x1[iptr8];
  
        UVW_y = a0*b0*c0*x1[iptr1+siz_bx] + a1*b0*c0*x2[iptr2+siz_bx]
               +a1*b1*c0*x2[iptr3+siz_bx] + a0*b1*c0*x1[iptr4+siz_bx]
               +a0*b0*c1*x1[iptr5+siz_bx] + a1*b0*c1*x2[iptr6+siz_bx]
               +a1*b1*c1*x2[iptr7+siz_bx] + a0*b1*c1*x1[iptr8+siz_bx];
  
        UVW_z = a0*b0*c0*x1[iptr1+2*siz_bx] + a1*b0*c0*x2[iptr2+2*siz_bx]
               +a1*b1*c0*x2[iptr3+2*siz_bx] + a0*b1*c0*x1[iptr4+2*siz_bx]
               +a0*b0*c1*x1[iptr5+2*siz_bx] + a1*b0*c1*x2[iptr6+2*siz_bx]
               +a1*b1*c1*x2[iptr7+2*siz_bx] + a0*b1*c1*x1[iptr8+2*siz_bx];
  
        iptr = k*siz_iz + j*siz_iy + i;  //(i,j,k) 
        x3d[iptr] = U_x + V_x + W_x - UV_x - UW_x - VW_x + UVW_x;
        y3d[iptr] = U_y + V_y + W_y - UV_y - UW_y - VW_y + UVW_y;
        z3d[iptr] = U_z + V_z + W_z - UV_z - UW_z - VW_z + UVW_z;
      }
    }
  }

  // assgn 6 bdry point to ghost
  // bdry x1
  if(mympi->neighid[0] == MPI_PROC_NULL)
  {
    for (int k=0; k<nz; k++)
    {
      for (int j=0; j<ny; j++)
      {
        iptr = k*siz_iz + j*siz_iy;  //(0,j,k) 
        gnk = gnk1 + k;
        gnj = gnj1 + j;
        iptr1 = gnk*total_ny + gnj;

        x3d[iptr] = x1[iptr1];
        y3d[iptr] = x1[iptr1+siz_bx];
        z3d[iptr] = x1[iptr1+2*siz_bx];
      }
    }
  }
  // bdry x2
  if(mympi->neighid[1] == MPI_PROC_NULL)
  {
    for (int k=0; k<nz; k++)
    {
      for (int j=0; j<ny; j++)
      {
        iptr = k*siz_iz + j*siz_iy + nx-1;  //(nx-1,j,k) 
        gnk = gnk1 + k;
        gnj = gnj1 + j;
        iptr1 = gnk*total_ny + gnj;

        x3d[iptr] = x2[iptr1];
        y3d[iptr] = x2[iptr1+siz_bx];
        z3d[iptr] = x2[iptr1+2*siz_bx];
      }
    }
  }
  // bdry y1
  if(mympi->neighid[2] == MPI_PROC_NULL)
  {
    for (int k=0; k<nz; k++)
    {
      for (int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + i;  //(i,0,k) 
        gnk = gnk1 + k; 
        gni = gni1 + i; 
        iptr1 = gnk*total_nx + gni;

        x3d[iptr] = y1[iptr1];
        y3d[iptr] = y1[iptr1+siz_by];
        z3d[iptr] = y1[iptr1+2*siz_by];
      }
    }
  }
  // bdry y2
  if(mympi->neighid[3] == MPI_PROC_NULL)
  {
    for (int k=0; k<nz; k++)
    {
      for (int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + (ny-1)*siz_iy + i;  //(i,ny-1,k) 
        gnk = gnk1 + k; 
        gni = gni1 + i; 
        iptr1 = gnk*total_nx + gni;

        x3d[iptr] = y2[iptr1];
        y3d[iptr] = y2[iptr1+siz_by];
        z3d[iptr] = y2[iptr1+2*siz_by];
      }
    }
  }
  // bdry z1
  if(mympi->neighid[4] == MPI_PROC_NULL)
  {
    for (int j=0; j<ny; j++)
    {
      for (int i=0; i<nx; i++)
      {
        iptr = j*siz_iy + i;  //(i,j,0) 
        gnj = gnj1 + j; 
        gni = gni1 + i; 
        iptr1 = gnj*total_nx + gni;

        x3d[iptr] = z1[iptr1];
        y3d[iptr] = z1[iptr1+siz_bz];
        z3d[iptr] = z1[iptr1+2*siz_bz];
      }
    }
  }
  // bdry z2
  if(mympi->neighid[5] == MPI_PROC_NULL)
  {
    for (int j=0; j<ny; j++)
    {
      for (int i=0; i<nx; i++)
      {
        iptr = (nz-1)*siz_iz + j*siz_iy + i;  //(i,j,nz-1) 
        gnj = gnj1 + j; 
        gni = gni1 + i; 
        iptr1 = gnj*total_nx + gni;

        x3d[iptr] = z2[iptr1];
        y3d[iptr] = z2[iptr1+siz_bz];
        z3d[iptr] = z2[iptr1+2*siz_bz];
      }
    }
  }
  
  return 0;
}

