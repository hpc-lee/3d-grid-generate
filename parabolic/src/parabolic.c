#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "parabolic.h"
#include "solver.h"
#include "lib_mem.h"

int 
para_gene(gd_t *gdcurv, mympi_t *mympi, bdry_t *bdry, par_t *par)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;

  float *x1 = bdry->x1;
  float *x2 = bdry->x2;
  float *y1 = bdry->y1;
  float *y2 = bdry->y2;

  MPI_Comm topocomm = mympi->topocomm;
  int myid = mympi->myid;
  int *neighid = mympi->neighid;

  int num_of_s_reqs = 2;
  int num_of_r_reqs = 2;

  float coef = par->coef;
  int o2i = par->o2i;

  bdry_effct_t *bdry_effct = (bdry_effct_t *) malloc(sizeof(bdry_effct_t)); 
  init_bdry_effct(bdry_effct,gdcurv);

  if(myid == 0 && o2i == 1)
  {
    fprintf(stdout,"parabolic method, march direction from inner bdry(k=0) to outer bdry(nz-1)\n");
    fprintf(stdout,"so if we want from outer bdry(nz-1) to inner bdry(k=0), must be flip\n");
    fprintf(stdout,"flip coord, trans outer bdry(nz-1) to inner bdry(0)\n");
  }
  if(o2i == 1)
  {
    flip_coord_z(gdcurv);
    flip_bdry_z(x1,x2,y1,y2,total_nx,total_ny,total_nz);
  }
  // malloc space for thomas method 
  // a,b,c,d_x,d_y,d_z,u_x,u_y,u_z. 9 vars space
  float *var_th = (float *)mem_calloc_1d_float((total_nx-2)*9, 0.0, "init");

  // for update
  float *coord = (float *) malloc(sizeof(float)*nx*ny*3);

  // calculate x1 x2 y1 y2 bdry arc length
  float *x1_len = (float *)mem_calloc_1d_float(total_nz*total_ny, 0.0, "init");
  float *x2_len = (float *)mem_calloc_1d_float(total_nz*total_ny, 0.0, "init");
  float *y1_len = (float *)mem_calloc_1d_float(total_nz*total_nx, 0.0, "init");
  float *y2_len = (float *)mem_calloc_1d_float(total_nz*total_nx, 0.0, "init");

  cal_bdry_arc_length(x1,x2,y1,y2,total_nx,total_ny,total_nz,x1_len,x2_len,y1_len,y2_len);

  for(int k=nk1; k<=nk2; k++)
  {
    // predict k, k+1 layer points
    predict_point(gdcurv,bdry_effct,mympi,k,o2i,coef,x1_len,x2_len,y1_len,y2_len);
    exchange_coord(gdcurv, mympi, k, num_of_s_reqs, num_of_r_reqs);
    // base predict points
    // update k layer points
    update_point(gdcurv,var_th,k,coord);
    exchange_coord(gdcurv, mympi, k, num_of_s_reqs, num_of_r_reqs);
    if(myid == 0)
    {
      fprintf(stdout,"number of layer is %d\n",k);
      fflush(stdout);
    }
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
  free(coord);

  return 0;
}

int 
predict_point(gd_t *gdcurv, bdry_effct_t *bdry_effct, mympi_t *mympi,
              int k, int o2i, float coef, float *x1_len,
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

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int gni1 = gdcurv->gni1;
  int gnj1 = gdcurv->gnj1;
  int gni, gnj;

  zt = (1.0*(k-1))/(nz-1-2);
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
      vn_x = vn_x/vn_len;
      vn_y = vn_y/vn_len;
      vn_z = vn_z/vn_len;
      
      // inner point i
      iptr1 = (k-1)*siz_iz + j*siz_iy + i; // (i,j,k-1)
      iptr2 = (nz-1)*siz_iz + j*siz_iy + i; // (i,j,nz-1)
      R_x = x3d[iptr1] - x3d[iptr2];
      R_y = y3d[iptr1] - y3d[iptr2];
      R_z = z3d[iptr1] - z3d[iptr2];
      R = sqrt(pow(R_x,2)+pow(R_y,2)+pow(R_z,2));

      // bdry point i=0
      gni = gni1 + i;
      gnj = gnj1 + j;
      iptr1 = (k-1)*total_ny + gnj;
      iptr2 = (nz-1)*total_ny + gnj;
      iptr3 = (k+1)*total_ny + gnj;
      iptr4 = k*total_ny + gnj;
      R_x1  = x1_len[iptr2] - x1_len[iptr1];
      r_x1  = x1_len[iptr3] - x1_len[iptr1];
      r_x11 = x1_len[iptr4] - x1_len[iptr1];

      c_x1  = r_x1/R_x1;
      c_x11 = r_x11/r_x1;

      // bdry point i=nx-1
      R_x2  = x2_len[iptr2] - x2_len[iptr1];
      r_x2  = x2_len[iptr3] - x2_len[iptr1];
      r_x22 = x2_len[iptr4] - x2_len[iptr1];

      c_x2  = r_x2/R_x2;
      c_x22 = r_x22/r_x2;

      // cal clustering factor c
      xi = (1.0*gni)/(total_nx-1);
      c_x =  (1-xi)*c_x1  + xi*c_x2;
      cc_x = (1-xi)*c_x11 + xi*c_x22;
      
      // bdry point j=0
      gni = gni1 + i;
      gnj = gnj1 + j;
      iptr1 = (k-1)*total_nx + gni;
      iptr2 = (nz-1)*total_nx + gni;
      iptr3 = (k+1)*total_nx + gni;
      iptr4 = k*total_nx + gni;
      R_y1  = y1_len[iptr2] - y1_len[iptr1];
      r_y1  = y1_len[iptr3] - y1_len[iptr1];
      r_y11 = y1_len[iptr4] - y1_len[iptr1];

      c_y1  = r_y1/R_y1;
      c_y11 = r_y11/r_y1;

      // bdry point j=ny-1
      R_y2  = y2_len[iptr2] - y2_len[iptr1];
      r_y2  = y2_len[iptr3] - y2_len[iptr1];
      r_y22 = y2_len[iptr4] - y2_len[iptr1];

      c_y2  = r_y2/R_y2;
      c_y22 = r_y22/r_y2;

      // cal clustering factor c
      et = (1.0*gnj)/(total_ny-1);
      c_y =  (1-et)*c_y1  + et*c_y2;
      cc_y = (1-et)*c_y11 + et*c_y22;

      c = 0.5*(c_x+c_y);
      cc = 0.5*(cc_x+cc_y);

      iptr1 = (k-1)*siz_iz + j*siz_iy + i; //(i,j,k-1)
      // x0,y0,z0 normal point
      x0 = x3d[iptr1] + vn_x*c*R;
      y0 = y3d[iptr1] + vn_y*c*R;
      z0 = z3d[iptr1] + vn_z*c*R;

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

  bdry_modify(gdcurv, bdry_effct, mympi, k);

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

    // i=nx-2 modify c,d_x,d_z
    iptr = k*siz_iz + j*siz_iy + (nx-1);
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

  return 0;
}

int
bdry_modify(gd_t *gdcurv, bdry_effct_t *bdry_effct, mympi_t *mympi, int k)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int gni1 =gdcurv->gni1; 
  int gnj1 =gdcurv->gnj1; 
  int gni, gnj;

  int *neighid = mympi->neighid;
  MPI_Comm topocomm = mympi->topocomm;

  // add bdry point effct
  // due to bdry maybe nonlinear,
  // bdry undulate
  size_t iptr1,iptr2,iptr3;

  float *dif_b1_x_k_loc  = bdry_effct->dif_b1_x_k_loc; 
  float *dif_b1_y_k_loc  = bdry_effct->dif_b1_y_k_loc; 
  float *dif_b1_z_k_loc  = bdry_effct->dif_b1_z_k_loc; 
  float *dif_b1_x_k1_loc = bdry_effct->dif_b1_x_k1_loc;
  float *dif_b1_y_k1_loc = bdry_effct->dif_b1_y_k1_loc;
  float *dif_b1_z_k1_loc = bdry_effct->dif_b1_z_k1_loc;

  float *dif_b1_x_k  = bdry_effct->dif_b1_x_k; 
  float *dif_b1_y_k  = bdry_effct->dif_b1_y_k; 
  float *dif_b1_z_k  = bdry_effct->dif_b1_z_k; 
  float *dif_b1_x_k1 = bdry_effct->dif_b1_x_k1;
  float *dif_b1_y_k1 = bdry_effct->dif_b1_y_k1;
  float *dif_b1_z_k1 = bdry_effct->dif_b1_z_k1;
  float b1_x1_k,  b1_y1_k,  b1_z1_k;
  float b1_x1_k1, b1_y1_k1, b1_z1_k1;
  float b1_x2_k,  b1_y2_k,  b1_z2_k;
  float b1_x2_k1, b1_y2_k1, b1_z2_k1;

  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++)
    {
      gnj = gnj1 + j;
      iptr1 = k*siz_iz + j*siz_iy + 0; //(0,j,k)
      iptr2 = (k-1)*siz_iz + j*siz_iy + 0; //(0,j,k-1)
      iptr3 = (k+1)*siz_iz + j*siz_iy + 0; //(0,j,k+1)

      b1_x1_k = x3d[iptr1]-x3d[iptr2];
      b1_y1_k = y3d[iptr1]-y3d[iptr2];
      b1_z1_k = z3d[iptr1]-z3d[iptr2];

      b1_x1_k1 = x3d[iptr3]-x3d[iptr1];
      b1_y1_k1 = y3d[iptr3]-y3d[iptr1];
      b1_z1_k1 = z3d[iptr3]-z3d[iptr1];

      iptr1 = k*siz_iz + j*siz_iy + 1; //(1,j,k)
      iptr2 = (k-1)*siz_iz + j*siz_iy + 1; //(1,j,k-1)
      iptr3 = (k+1)*siz_iz + j*siz_iy + 1; //(1,j,k+1)

      b1_x2_k = x3d[iptr1]-x3d[iptr2];
      b1_y2_k = y3d[iptr1]-y3d[iptr2];
      b1_z2_k = z3d[iptr1]-z3d[iptr2];

      b1_x2_k1 = x3d[iptr3]-x3d[iptr1];
      b1_y2_k1 = y3d[iptr3]-y3d[iptr1];
      b1_z2_k1 = z3d[iptr3]-z3d[iptr1];

      dif_b1_x_k_loc[gnj]  = b1_x1_k -b1_x2_k;
      dif_b1_y_k_loc[gnj]  = b1_y1_k -b1_y2_k;
      dif_b1_z_k_loc[gnj]  = b1_z1_k -b1_z2_k;
      dif_b1_x_k1_loc[gnj] = b1_x1_k1-b1_x2_k1;
      dif_b1_y_k1_loc[gnj] = b1_y1_k1-b1_y2_k1;
      dif_b1_z_k1_loc[gnj] = b1_z1_k1-b1_z2_k1;
    }
  }

  float *dif_b2_x_k_loc  = bdry_effct->dif_b2_x_k_loc; 
  float *dif_b2_y_k_loc  = bdry_effct->dif_b2_y_k_loc; 
  float *dif_b2_z_k_loc  = bdry_effct->dif_b2_z_k_loc; 
  float *dif_b2_x_k1_loc = bdry_effct->dif_b2_x_k1_loc;
  float *dif_b2_y_k1_loc = bdry_effct->dif_b2_y_k1_loc;
  float *dif_b2_z_k1_loc = bdry_effct->dif_b2_z_k1_loc;

  float *dif_b2_x_k  = bdry_effct->dif_b2_x_k; 
  float *dif_b2_y_k  = bdry_effct->dif_b2_y_k; 
  float *dif_b2_z_k  = bdry_effct->dif_b2_z_k; 
  float *dif_b2_x_k1 = bdry_effct->dif_b2_x_k1;
  float *dif_b2_y_k1 = bdry_effct->dif_b2_y_k1;
  float *dif_b2_z_k1 = bdry_effct->dif_b2_z_k1;
  float b2_x1_k,  b2_y1_k,  b2_z1_k;
  float b2_x1_k1, b2_y1_k1, b2_z1_k1;
  float b2_x2_k,  b2_y2_k,  b2_z2_k;
  float b2_x2_k1, b2_y2_k1, b2_z2_k1;

  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++)
    {
      gnj = gnj1 + j;
      iptr1 = k*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k)
      iptr2 = (k-1)*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k-1)
      iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k+1)

      b2_x1_k  = x3d[iptr1]-x3d[iptr2];
      b2_y1_k  = y3d[iptr1]-y3d[iptr2];
      b2_z1_k  = z3d[iptr1]-z3d[iptr2];

      b2_x1_k1 = x3d[iptr3]-x3d[iptr1];
      b2_y1_k1 = y3d[iptr3]-y3d[iptr1];
      b2_z1_k1 = z3d[iptr3]-z3d[iptr1];

      iptr1 = k*siz_iz + j*siz_iy + nx-2; //(nx-2,j,k)
      iptr2 = (k-1)*siz_iz + j*siz_iy + nx-2; //(nx-2,j,k-1)
      iptr3 = (k+1)*siz_iz + j*siz_iy + nx-2; //(nx-2,j,k+1)

      b2_x2_k  = x3d[iptr1]-x3d[iptr2];
      b2_y2_k  = y3d[iptr1]-y3d[iptr2];
      b2_z2_k  = z3d[iptr1]-z3d[iptr2];

      b2_x2_k1 = x3d[iptr3]-x3d[iptr1];
      b2_y2_k1 = y3d[iptr3]-y3d[iptr1];
      b2_z2_k1 = z3d[iptr3]-z3d[iptr1];

      dif_b2_x_k_loc[gnj]  = b2_x1_k -b2_x2_k;
      dif_b2_y_k_loc[gnj]  = b2_y1_k -b2_y2_k;
      dif_b2_z_k_loc[gnj]  = b2_z1_k -b2_z2_k;
      dif_b2_x_k1_loc[gnj] = b2_x1_k1-b2_x2_k1;
      dif_b2_y_k1_loc[gnj] = b2_y1_k1-b2_y2_k1;
      dif_b2_z_k1_loc[gnj] = b2_z1_k1-b2_z2_k1;
    }
  }

  float *dif_b3_x_k_loc  = bdry_effct->dif_b3_x_k_loc; 
  float *dif_b3_y_k_loc  = bdry_effct->dif_b3_y_k_loc; 
  float *dif_b3_z_k_loc  = bdry_effct->dif_b3_z_k_loc; 
  float *dif_b3_x_k1_loc = bdry_effct->dif_b3_x_k1_loc;
  float *dif_b3_y_k1_loc = bdry_effct->dif_b3_y_k1_loc;
  float *dif_b3_z_k1_loc = bdry_effct->dif_b3_z_k1_loc;

  float *dif_b3_x_k  = bdry_effct->dif_b3_x_k; 
  float *dif_b3_y_k  = bdry_effct->dif_b3_y_k; 
  float *dif_b3_z_k  = bdry_effct->dif_b3_z_k; 
  float *dif_b3_x_k1 = bdry_effct->dif_b3_x_k1;
  float *dif_b3_y_k1 = bdry_effct->dif_b3_y_k1;
  float *dif_b3_z_k1 = bdry_effct->dif_b3_z_k1;
  float b3_x1_k,  b3_y1_k,  b3_z1_k;
  float b3_x1_k1, b3_y1_k1, b3_z1_k1;
  float b3_x2_k,  b3_y2_k,  b3_z2_k;
  float b3_x2_k1, b3_y2_k1, b3_z2_k1;

  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      gni = gni1 + i;
      iptr1 = k*siz_iz + 0*siz_iy + i; //(i,0,k)
      iptr2 = (k-1)*siz_iz + 0*siz_iy + i; //(i,0,k-1)
      iptr3 = (k+1)*siz_iz + 0*siz_iy + i; //(i,0,k+1)

      b3_x1_k  = x3d[iptr1]-x3d[iptr2];
      b3_y1_k  = y3d[iptr1]-y3d[iptr2];
      b3_z1_k  = z3d[iptr1]-z3d[iptr2];

      b3_x1_k1 = x3d[iptr3]-x3d[iptr1];
      b3_y1_k1 = y3d[iptr3]-y3d[iptr1];
      b3_z1_k1 = z3d[iptr3]-z3d[iptr1];

      iptr1 = k*siz_iz + 1*siz_iy + i; //(i,1,k)
      iptr2 = (k-1)*siz_iz + 1*siz_iy + i; //(i,1,k-1)
      iptr3 = (k+1)*siz_iz + 1*siz_iy + i; //(i,1,k+1)

      b3_x2_k = x3d[iptr1]-x3d[iptr2];
      b3_y2_k = y3d[iptr1]-y3d[iptr2];
      b3_z2_k = z3d[iptr1]-z3d[iptr2];

      b3_x2_k1 = x3d[iptr3]-x3d[iptr1];
      b3_y2_k1 = y3d[iptr3]-y3d[iptr1];
      b3_z2_k1 = z3d[iptr3]-z3d[iptr1];

      dif_b3_x_k_loc[gni]  = b3_x1_k -b3_x2_k;
      dif_b3_y_k_loc[gni]  = b3_y1_k -b3_y2_k;
      dif_b3_z_k_loc[gni]  = b3_z1_k -b3_z2_k;
      dif_b3_x_k1_loc[gni] = b3_x1_k1-b3_x2_k1;
      dif_b3_y_k1_loc[gni] = b3_y1_k1-b3_y2_k1;
      dif_b3_z_k1_loc[gni] = b3_z1_k1-b3_z2_k1;
    }
  }

  float *dif_b4_x_k_loc  = bdry_effct->dif_b4_x_k_loc; 
  float *dif_b4_y_k_loc  = bdry_effct->dif_b4_y_k_loc; 
  float *dif_b4_z_k_loc  = bdry_effct->dif_b4_z_k_loc; 
  float *dif_b4_x_k1_loc = bdry_effct->dif_b4_x_k1_loc;
  float *dif_b4_y_k1_loc = bdry_effct->dif_b4_y_k1_loc;
  float *dif_b4_z_k1_loc = bdry_effct->dif_b4_z_k1_loc;

  float *dif_b4_x_k  = bdry_effct->dif_b4_x_k; 
  float *dif_b4_y_k  = bdry_effct->dif_b4_y_k; 
  float *dif_b4_z_k  = bdry_effct->dif_b4_z_k; 
  float *dif_b4_x_k1 = bdry_effct->dif_b4_x_k1;
  float *dif_b4_y_k1 = bdry_effct->dif_b4_y_k1;
  float *dif_b4_z_k1 = bdry_effct->dif_b4_z_k1;
  float b4_x1_k,  b4_y1_k,  b4_z1_k;
  float b4_x1_k1, b4_y1_k1, b4_z1_k1;
  float b4_x2_k,  b4_y2_k,  b4_z2_k;
  float b4_x2_k1, b4_y2_k1, b4_z2_k1;

  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      gni = gni1 + i;
      iptr1 = k*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k)
      iptr2 = (k-1)*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k-1)
      iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k+1)

      b4_x1_k = x3d[iptr1]-x3d[iptr2];
      b4_y1_k = y3d[iptr1]-y3d[iptr2];
      b4_z1_k = z3d[iptr1]-z3d[iptr2];

      b4_x1_k1 = x3d[iptr3]-x3d[iptr1];
      b4_y1_k1 = y3d[iptr3]-y3d[iptr1];
      b4_z1_k1 = z3d[iptr3]-z3d[iptr1];

      iptr1 = k*siz_iz + (ny-2)*siz_iy + i; //(i,ny-2,k)
      iptr2 = (k-1)*siz_iz + (ny-2)*siz_iy + i; //(i,ny-2,k-1)
      iptr3 = (k+1)*siz_iz + (ny-2)*siz_iy + i; //(i,ny-2,k+1)

      b4_x2_k = x3d[iptr1]-x3d[iptr2];
      b4_y2_k = y3d[iptr1]-y3d[iptr2];
      b4_z2_k = z3d[iptr1]-z3d[iptr2];

      b4_x2_k1 = x3d[iptr3]-x3d[iptr1];
      b4_y2_k1 = y3d[iptr3]-y3d[iptr1];
      b4_z2_k1 = z3d[iptr3]-z3d[iptr1];

      dif_b4_x_k_loc[gni]  = b4_x1_k -b4_x2_k;
      dif_b4_y_k_loc[gni]  = b4_y1_k -b4_y2_k;
      dif_b4_z_k_loc[gni]  = b4_z1_k -b4_z2_k;
      dif_b4_x_k1_loc[gni] = b4_x1_k1-b4_x2_k1;
      dif_b4_y_k1_loc[gni] = b4_y1_k1-b4_y2_k1;
      dif_b4_z_k1_loc[gni] = b4_z1_k1-b4_z2_k1;
    }
  }

  MPI_Barrier(topocomm);

  MPI_Allreduce(dif_b1_x_k_loc,  dif_b1_x_k,  total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b1_y_k_loc,  dif_b1_y_k,  total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b1_z_k_loc,  dif_b1_z_k,  total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b1_x_k1_loc, dif_b1_x_k1, total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b1_y_k1_loc, dif_b1_y_k1, total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b1_z_k1_loc, dif_b1_z_k1, total_ny, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(dif_b2_x_k_loc,  dif_b2_x_k,  total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b2_y_k_loc,  dif_b2_y_k,  total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b2_z_k_loc,  dif_b2_z_k,  total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b2_x_k1_loc, dif_b2_x_k1, total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b2_y_k1_loc, dif_b2_y_k1, total_ny, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b2_z_k1_loc, dif_b2_z_k1, total_ny, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(dif_b3_x_k_loc,  dif_b3_x_k,  total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b3_y_k_loc,  dif_b3_y_k,  total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b3_z_k_loc,  dif_b3_z_k,  total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b3_x_k1_loc, dif_b3_x_k1, total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b3_y_k1_loc, dif_b3_y_k1, total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b3_z_k1_loc, dif_b3_z_k1, total_nx, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(dif_b4_x_k_loc,  dif_b4_x_k,  total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b4_y_k_loc,  dif_b4_y_k,  total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b4_z_k_loc,  dif_b4_z_k,  total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b4_x_k1_loc, dif_b4_x_k1, total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b4_y_k1_loc, dif_b4_y_k1, total_nx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(dif_b4_z_k1_loc, dif_b4_z_k1, total_nx, MPI_FLOAT, MPI_SUM, topocomm);

  float xi, et;
  float bdry_x_k, bdry_y_k, bdry_z_k;
  float bdry_x_k1, bdry_y_k1, bdry_z_k1;
  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      gni = gni1 + i;
      gnj = gnj1 + j;
      xi = (1.0*gni)/(total_nx-1);
      bdry_x_k = (1-xi)*dif_b1_x_k[gnj] + xi*dif_b2_x_k[gnj];
      bdry_y_k = (1-xi)*dif_b1_y_k[gnj] + xi*dif_b2_y_k[gnj];
      bdry_z_k = (1-xi)*dif_b1_z_k[gnj] + xi*dif_b2_z_k[gnj];
      iptr1 = k*siz_iz + j*siz_iy + i; 
      x3d[iptr1] = x3d[iptr1] + 1.32*bdry_x_k;
      y3d[iptr1] = y3d[iptr1] + 1.32*bdry_y_k;
      z3d[iptr1] = z3d[iptr1] + 1.32*bdry_z_k;

      bdry_x_k1 = (1-xi)*dif_b1_x_k1[gnj] + xi*dif_b2_x_k1[gnj];
      bdry_y_k1 = (1-xi)*dif_b1_y_k1[gnj] + xi*dif_b2_y_k1[gnj];
      bdry_z_k1 = (1-xi)*dif_b1_z_k1[gnj] + xi*dif_b2_z_k1[gnj];
      iptr2 = (k+1)*siz_iz + j*siz_iy + i; 
      x3d[iptr2] = x3d[iptr2] + 1.32*bdry_x_k1;
      y3d[iptr2] = y3d[iptr2] + 1.32*bdry_y_k1;
      z3d[iptr2] = z3d[iptr2] + 1.32*bdry_z_k1;
    }
  }

  for(int j=1; j<ny-1; j++) {
    for(int i=1; i<nx-1; i++)
    {
      gni = gni1 + i;
      gnj = gnj1 + j;
      et = (1.0*gnj)/(total_ny-1);
      bdry_x_k = (1-et)*dif_b3_x_k[gni] + et*dif_b4_x_k[gni];
      bdry_y_k = (1-et)*dif_b3_y_k[gni] + et*dif_b4_y_k[gni];
      bdry_z_k = (1-et)*dif_b3_z_k[gni] + et*dif_b4_z_k[gni];
      iptr1 = k*siz_iz + j*siz_iy + i;
      x3d[iptr1] = x3d[iptr1] + 1.30*bdry_x_k;
      y3d[iptr1] = y3d[iptr1] + 1.30*bdry_y_k;
      z3d[iptr1] = z3d[iptr1] + 1.30*bdry_z_k;

      bdry_x_k1 = (1-et)*dif_b3_x_k1[gni] + et*dif_b4_x_k1[gni];
      bdry_y_k1 = (1-et)*dif_b3_y_k1[gni] + et*dif_b4_y_k1[gni];
      bdry_z_k1 = (1-et)*dif_b3_z_k1[gni] + et*dif_b4_z_k1[gni];
      iptr2 = (k+1)*siz_iz + j*siz_iy + i;
      x3d[iptr2] = x3d[iptr2] + 1.57*bdry_x_k1;
      y3d[iptr2] = y3d[iptr2] + 1.57*bdry_y_k1;
      z3d[iptr2] = z3d[iptr2] + 1.57*bdry_z_k1;
    }
  }

  return 0;
}

int 
init_bdry_effct(bdry_effct_t *bdry_effct, gd_t *gdcurv)
{
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;

  // x1 bdry1
  bdry_effct->dif_b1_x_k_loc  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_y_k_loc  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_z_k_loc  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_x_k1_loc = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_y_k1_loc = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_z_k1_loc = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");

  bdry_effct->dif_b1_x_k  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_y_k  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_z_k  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_x_k1 = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_y_k1 = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b1_z_k1 = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");

  // x2 bdry2
  bdry_effct->dif_b2_x_k_loc  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_y_k_loc  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_z_k_loc  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_x_k1_loc = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_y_k1_loc = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_z_k1_loc = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");

  bdry_effct->dif_b2_x_k  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_y_k  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_z_k  = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_x_k1 = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_y_k1 = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");
  bdry_effct->dif_b2_z_k1 = (float *)mem_calloc_1d_float(total_ny, 0.0, "init");

  // y1 bdry3
  bdry_effct->dif_b3_x_k_loc  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_y_k_loc  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_z_k_loc  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_x_k1_loc = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_y_k1_loc = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_z_k1_loc = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");

  bdry_effct->dif_b3_x_k  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_y_k  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_z_k  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_x_k1 = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_y_k1 = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b3_z_k1 = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");

  // y2 bdry4
  bdry_effct->dif_b4_x_k_loc  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_y_k_loc  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_z_k_loc  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_x_k1_loc = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_y_k1_loc = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_z_k1_loc = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");

  bdry_effct->dif_b4_x_k  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_y_k  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_z_k  = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_x_k1 = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_y_k1 = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");
  bdry_effct->dif_b4_z_k1 = (float *)mem_calloc_1d_float(total_nx, 0.0, "init");

  return 0;
}

int
flip_bdry_z(float *x1, float *x2, float *y1, float *y2,
            int total_nx, int total_ny, int total_nz)
{
  int size_bx = total_ny*total_nz;
  int size_by = total_nx*total_nz;
  float *tmp_x = (float *) malloc(3*size_bx*sizeof(float));
  float *tmp_y = (float *) malloc(3*size_by*sizeof(float));
  size_t iptr, iptr1;
  
  // x1
  // copy bdry
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      tmp_x[iptr]           = x1[iptr];
      tmp_x[iptr+size_bx]   = x1[iptr+size_bx];
      tmp_x[iptr+2*size_bx] = x1[iptr+2*size_bx];
    }
  }
  // flip bdry
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr  = k*total_ny + j;
      iptr1 = (total_nz-k-1)*total_ny + j;
      x1[iptr]           = tmp_x[iptr1];
      x1[iptr+size_bx]   = tmp_x[iptr1+size_bx];
      x1[iptr+2*size_bx] = tmp_x[iptr1+2*size_bx];
    }
  }

  // x2
  // copy bdry
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      tmp_x[iptr]           = x2[iptr];
      tmp_x[iptr+size_bx]   = x2[iptr+size_bx];
      tmp_x[iptr+2*size_bx] = x2[iptr+2*size_bx];
    }
  }
  // flip bdry
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr  = k*total_ny + j;
      iptr1 = (total_nz-k-1)*total_ny + j;
      x2[iptr]           = tmp_x[iptr1];
      x2[iptr+size_bx]   = tmp_x[iptr1+size_bx];
      x2[iptr+2*size_bx] = tmp_x[iptr1+2*size_bx];
    }
  }

  // y1
  // copy bdry
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      tmp_y[iptr]           = y1[iptr];
      tmp_y[iptr+size_by]   = y1[iptr+size_by];
      tmp_y[iptr+2*size_by] = y1[iptr+2*size_by];
    }
  }
  // flip bdry
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr  = k*total_nx + i;
      iptr1 = (total_nz-k-1)*total_nx + i;
      y1[iptr]           = tmp_y[iptr1];
      y1[iptr+size_by]   = tmp_y[iptr1+size_by];
      y1[iptr+2*size_by] = tmp_y[iptr1+2*size_by];
    }
  }

  // y2
  // copy bdry
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      tmp_y[iptr]           = y2[iptr];
      tmp_y[iptr+size_by]   = y2[iptr+size_by];
      tmp_y[iptr+2*size_by] = y2[iptr+2*size_by];
    }
  }
  // flip bdry
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr  = k*total_nx + i;
      iptr1 = (total_nz-k-1)*total_nx + i;
      y2[iptr]           = tmp_y[iptr1];
      y2[iptr+size_by]   = tmp_y[iptr1+size_by];
      y2[iptr+2*size_by] = tmp_y[iptr1+2*size_by];
    }
  }

  free(tmp_x);
  free(tmp_y);

  return 0;
}

int
cal_bdry_arc_length(float *x1, float *x2, float *y1, float *y2, 
                    int total_nx, int total_ny, int total_nz, 
                    float *x1_len, float *x2_len, float *y1_len, float *y2_len)
{
  size_t iptr1, iptr2, iptr3, iptr4;
  float x_len, y_len, z_len, dh_len;

  size_t size_bx = total_ny*total_nz;
  size_t size_by = total_nx*total_nz;
  for(int k=1; k<total_nz; k++) {
    for(int j=1; j<total_ny-1; j++)
    {
      iptr1 = k*total_ny + j;
      iptr2 = (k-1)*total_ny + j;

      x_len = x1[iptr1+0*size_bx] - x1[iptr2+0*size_bx];
      y_len = x1[iptr1+1*size_bx] - x1[iptr2+1*size_bx];
      z_len = x1[iptr1+2*size_bx] - x1[iptr2+2*size_bx];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      x1_len[iptr1] = x1_len[iptr2] + dh_len;

      x_len = x2[iptr1+0*size_bx] - x2[iptr2+0*size_bx];
      y_len = x2[iptr1+1*size_bx] - x2[iptr2+1*size_bx];
      z_len = x2[iptr1+2*size_bx] - x2[iptr2+2*size_bx];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      x2_len[iptr1] = x2_len[iptr2] + dh_len;
    }
  }

  for(int k=1; k<total_nz; k++) {
    for(int i=1; i<total_nx-1; i++)
    {
      iptr1 = k*total_nx + i;
      iptr2 = (k-1)*total_nx + i;

      x_len = y1[iptr1+0*size_by] - y1[iptr2+0*size_by];
      y_len = y1[iptr1+1*size_by] - y1[iptr2+1*size_by];
      z_len = y1[iptr1+2*size_by] - y1[iptr2+2*size_by];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      y1_len[iptr1] = y1_len[iptr2] + dh_len;

      x_len = y2[iptr1+0*size_by] - y2[iptr2+0*size_by];
      y_len = y2[iptr1+1*size_by] - y2[iptr2+1*size_by];
      z_len = y2[iptr1+2*size_by] - y2[iptr2+2*size_by];
      dh_len = sqrt(pow(x_len,2) + pow(y_len,2) + pow(z_len,2));
      y2_len[iptr1] = y2_len[iptr2] + dh_len;
    }
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

