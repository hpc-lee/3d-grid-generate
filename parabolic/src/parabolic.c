#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "parabolic.h"
#include "solver.h"
#include "lib_mem.h"

int 
para_gene(gd_t *gdcurv, mympi_t *mympi, par_t *par)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;
  float *step = gdcurv->step;

  MPI_Comm topocomm = mympi->topocomm;
  int myid = mympi->myid;
  int *neighid = mympi->neighid;

  int num_of_s_reqs = 2;
  int num_of_r_reqs = 2;

  float coef = par->coef;
  int t2b = par->t2b;

  if(myid == 0 && t2b == 1)
  {
    fprintf(stdout,"march direction from top bdry(k=nz-1) to bottom bdry(k=1)\n");
    fprintf(stdout,"so we flip bdry, set top bdry to bottom bdry\n");
  }
  if(t2b == 1)
  {
    flip_coord_z(gdcurv);
    flip_step_z(step,nz);
  }
  // malloc space for thomas method 
  // a,b,c,d_x,d_y,d_z,u_x,u_y,u_z. 9 vars space
  float *var_th = (float *)mem_calloc_1d_float((total_nx-2)*9, 0.0, "init");
  // for update
  float *coord = (float *) malloc(sizeof(float)*nx*ny*3);

  // calculate step_length
  float *step_len = (float *)mem_calloc_1d_float(nz, 0.0, "init");
  for(int k=1; k<nz; k++)
  {
    step_len[k] = step_len[k-1] + step[k-1];
  }

  for(int k=1; k<nz-1; k++)
  {
    // predict k, k+1 layer points
    predict_point(gdcurv,mympi,k,t2b,coef,step_len);
    exchange_coord(gdcurv,mympi,k,num_of_s_reqs,num_of_r_reqs);
    // base predict points
    // update k layer points
    update_point(gdcurv,var_th,k,coord);
    exchange_coord(gdcurv, mympi, k, num_of_s_reqs, num_of_r_reqs);
    if(myid == 0)
    {
      fprintf(stdout,"number of layer is %d\n",k+1);
      fflush(stdout);
    }
  }

  if(t2b == 1)
  {
    // flip to return.
    flip_coord_z(gdcurv);
  }

  free(var_th);
  free(coord);

  return 0;
}

int 
predict_point(gd_t *gdcurv, mympi_t *mympi, int k, int t2b,
              float coef, float *step_len)
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
  float R_x,R_y,R_z,R;
  float R1,r1,r2;
  float c1,c2;
  float x0,y0,z0,xs,ys,zs;
  float x_xi, y_xi, z_xi;
  float x_et, y_et, z_et;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;

  zt = (1.0*k)/(nz-1);
  cs = exp(-coef*zt);

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
      if(t2b == 1)
      {
        vn_x = -vn_x;
        vn_y = -vn_y;
        vn_z = -vn_z;
      } 
      vn_len = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
      vn_x = vn_x/vn_len;
      vn_y = vn_y/vn_len;
      vn_z = vn_z/vn_len;
      
      // inner point i
      iptr1 = (k-1)*siz_iz + j*siz_iy + i; // (i,j,k-1)
      iptr2 = (nz-1)*siz_iz + j*siz_iy + i; // (i,j,nz-1)
      R_x = x3d[iptr2] - x3d[iptr1];
      R_y = y3d[iptr2] - y3d[iptr1];
      R_z = z3d[iptr2] - z3d[iptr1];
      R = sqrt(pow(R_x,2)+pow(R_y,2)+pow(R_z,2));

      // cal clustering factor c
      R1 = step_len[nz-1] - step_len[k-1];
      r1 = step_len[k+1] - step_len[k-1];
      r2 = step_len[k] - step_len[k-1];
      c1 = r1/R1;
      c2 = r2/r1;
      
      // x0,y0,z0 normal point
      x0 = x3d[iptr1] + vn_x*c1*R;
      y0 = y3d[iptr1] + vn_y*c1*R;
      z0 = z3d[iptr1] + vn_z*c1*R;

      // xs,ys,zs linear distance point
      xs = x3d[iptr1] + c1*R_x;
      ys = y3d[iptr1] + c1*R_y;
      zs = z3d[iptr1] + c1*R_z;

      if(k<nz-2)
      {
        iptr3 = (k+1)*siz_iz + j*siz_iy + i;  //(i,j,k+1)
        x3d[iptr3] = cs*x0 + (1-cs)*xs;
        y3d[iptr3] = cs*y0 + (1-cs)*ys;
        z3d[iptr3] = cs*z0 + (1-cs)*zs;

        iptr4 = k*siz_iz + j*siz_iy + i;  //(i,j,k)
        x3d[iptr4] = x3d[iptr1] + c2*(x3d[iptr3]-x3d[iptr1]);
        y3d[iptr4] = y3d[iptr1] + c2*(y3d[iptr3]-y3d[iptr1]);
        z3d[iptr4] = z3d[iptr1] + c2*(z3d[iptr3]-z3d[iptr1]);
      }
      if(k==nz-2)
      {
        //bottom boundary is fixed
        iptr3 = (k+1)*siz_iz + j*siz_iy + i;  //(i,j,k+1)
        iptr4 = k*siz_iz + j*siz_iy + i;  //(i,j,k)
        x3d[iptr4] = x3d[iptr1] + c2*(x3d[iptr3]-x3d[iptr1]);
        y3d[iptr4] = y3d[iptr1] + c2*(y3d[iptr3]-y3d[iptr1]);
        z3d[iptr4] = z3d[iptr1] + c2*(z3d[iptr3]-z3d[iptr1]);
      }
    }
  }
  // geometric symmetry bdry
  for(int j=0; j<ny; j++) {
    iptr1 = (k+1)*siz_iz + j*siz_iy + 0;
    iptr2 = (k+1)*siz_iz + j*siz_iy + 1;
    iptr3 = (k+1)*siz_iz + j*siz_iy + 2;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];

    iptr1 = k*siz_iz + j*siz_iy + 0;
    iptr2 = k*siz_iz + j*siz_iy + 1;
    iptr3 = k*siz_iz + j*siz_iy + 2;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
    
    iptr1 = (k+1)*siz_iz + j*siz_iy + nx-1;
    iptr2 = (k+1)*siz_iz + j*siz_iy + nx-2;
    iptr3 = (k+1)*siz_iz + j*siz_iy + nx-3;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];

    iptr1 = k*siz_iz + j*siz_iy + nx-1;
    iptr2 = k*siz_iz + j*siz_iy + nx-2;
    iptr3 = k*siz_iz + j*siz_iy + nx-3;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
  }
  for(int i=0; i<nx; i++) {
    iptr1 = (k+1)*siz_iz + 0*siz_iy + i;
    iptr2 = (k+1)*siz_iz + 1*siz_iy + i;
    iptr3 = (k+1)*siz_iz + 2*siz_iy + i;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];

    iptr1 = k*siz_iz + 0*siz_iy + i;
    iptr2 = k*siz_iz + 1*siz_iy + i;
    iptr3 = k*siz_iz + 2*siz_iy + i;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
    
    iptr1 = (k+1)*siz_iz + (ny-1)*siz_iy + i;
    iptr2 = (k+1)*siz_iz + (ny-2)*siz_iy + i;
    iptr3 = (k+1)*siz_iz + (ny-3)*siz_iy + i;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];

    iptr1 = k*siz_iz + (ny-1)*siz_iy + i;
    iptr2 = k*siz_iz + (ny-2)*siz_iy + i;
    iptr3 = k*siz_iz + (ny-3)*siz_iy + i;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
  }

  return 0;
}

int
update_point(gd_t *gdcurv, float *var_th, int k, float *coord)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int siz_vec = nx-2;  // ni
  float *a = var_th;
  float *b = var_th + 1*siz_vec;
  float *c = var_th + 2*siz_vec;
  float *d_x = var_th + 3*siz_vec;
  float *d_y = var_th + 4*siz_vec;
  float *d_z = var_th + 5*siz_vec;
  float *u_x = var_th + 6*siz_vec;
  float *u_y = var_th + 7*siz_vec;
  float *u_z = var_th + 8*siz_vec;

  float *coord_x = coord;
  float *coord_y = coord + 1*nx*ny;
  float *coord_z = coord + 2*nx*ny;

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

  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = k*siz_iz + j*siz_iy + (i+1);  // (i+1,j,k)
      iptr2 = k*siz_iz + j*siz_iy + (i-1);  // (i-1,j,k)
      x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
      y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
      z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);

      iptr3 = k*siz_iz + (j+1)*siz_iy + i;  // (i,j+1,k)
      iptr4 = k*siz_iz + (j-1)*siz_iy + i;  // (i,j-1,k)
      x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
      y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
      z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);

      iptr5 = (k+1)*siz_iz + j*siz_iy + i;  // (i,j,k+1)
      iptr6 = (k-1)*siz_iz + j*siz_iy + i;  // (i,j,k-1)
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

      iptr1 = k*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,k)
      iptr2 = k*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,k)
      iptr3 = (k+1)*siz_iz + j*siz_iy + i;   //(i,j,k+1)
      iptr4 = (k-1)*siz_iz + j*siz_iy + i;   //(i,j,k-1)

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

    // i=nx-2 modify c,d_x,d_z
    iptr = k*siz_iz + j*siz_iy + (nx-1); // (nx-1,j,k)
    d_x[nx-3] = d_x[nx-3] - c[nx-3]*x3d[iptr];
    d_y[nx-3] = d_y[nx-3] - c[nx-3]*y3d[iptr];
    d_z[nx-3] = d_z[nx-3] - c[nx-3]*z3d[iptr];
    
    // cal coords and update
    thomas(siz_vec,a,b,c,d_x,d_y,d_z,u_x,u_y,u_z);

    for(int i=1; i<nx-1; i++)
    {
      iptr1 = j*siz_iy + i;
      coord_x[iptr1] = u_x[i-1];
      coord_y[iptr1] = u_y[i-1];
      coord_z[iptr1] = u_z[i-1];
    }
  }

  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr = k*siz_iz + j*siz_iy + i;
      iptr1 = j*siz_iy + i;
      x3d[iptr] = coord_x[iptr1];
      y3d[iptr] = coord_y[iptr1];
      z3d[iptr] = coord_z[iptr1];
    }
  }

  // geometric symmetry bdry
  for(int j=0; j<ny; j++) {
    iptr1 = k*siz_iz + j*siz_iy + 0;
    iptr2 = k*siz_iz + j*siz_iy + 1;
    iptr3 = k*siz_iz + j*siz_iy + 2;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
    
    iptr1 = k*siz_iz + j*siz_iy + nx-1;
    iptr2 = k*siz_iz + j*siz_iy + nx-2;
    iptr3 = k*siz_iz + j*siz_iy + nx-3;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
  }

  for(int i=0; i<nx; i++) {
    iptr1 = k*siz_iz + 0*siz_iy + i;
    iptr2 = k*siz_iz + 1*siz_iy + i;
    iptr3 = k*siz_iz + 2*siz_iy + i;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
    
    iptr1 = k*siz_iz + (ny-1)*siz_iy + i;
    iptr2 = k*siz_iz + (ny-2)*siz_iy + i;
    iptr3 = k*siz_iz + (ny-3)*siz_iy + i;
    x3d[iptr1] = 2*x3d[iptr2] - x3d[iptr3];
    y3d[iptr1] = 2*y3d[iptr2] - y3d[iptr3];
    z3d[iptr1] = 2*z3d[iptr2] - z3d[iptr3];
  }

  return 0;
}

int
exchange_coord(gd_t *gdcurv, mympi_t *mympi, int k, int num_of_s_reqs, int num_of_r_reqs)
{
  MPI_Startall(num_of_r_reqs, mympi->r_reqs);
  grid_pack_mesg(mympi,gdcurv,k);
  
  MPI_Startall(num_of_s_reqs, mympi->s_reqs);

  MPI_Waitall(num_of_s_reqs, mympi->s_reqs,MPI_STATUS_IGNORE);
  MPI_Waitall(num_of_r_reqs, mympi->r_reqs,MPI_STATUS_IGNORE);
  grid_unpack_mesg(mympi,gdcurv,k);

  return 0;
}

int
flip_step_z(float *step, int nz)
{
  float *step_tmp = (float *)malloc((nz-1)*sizeof(float));
  for(int k=0; k<nz-1; k++)
  {
    step_tmp[k] = step[k]; 
  }
  for(int k=0; k<nz-1; k++)
  {
    step[nz-2-k] = step_tmp[k]; 
  }

  free(step_tmp);

  return 0;
}

