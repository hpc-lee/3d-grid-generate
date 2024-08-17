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
hyper_gene(gd_t *gdcurv, par_t *par)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t size = (nx-2)*(ny-2);  // not include bdry 2 points

  float coef = par->coef;
  int o2i = par->o2i;
  int bdry_x_itype = par->bdry_x_itype;
  float epsilon_x = par->epsilon_x;
  int bdry_y_itype = par->bdry_y_itype;
  float epsilon_y = par->epsilon_y;

  float *coef_e_xi = NULL;
  coef_e_xi = (float *)mem_calloc_1d_float(size, 0.0, "init");
  float *coef_e_et = NULL;
  coef_e_et = (float *)mem_calloc_1d_float(size, 0.0, "init");
  float *volume = NULL;
  volume = (float *)mem_calloc_1d_float(nx*ny*2, 0.0, "init");
  // malloc space for thomas_block method 
  // each element is 3x3 matrix
  // xi
  double *a_xi = NULL;
  a_xi = (double *)mem_calloc_1d_double(size*3*3, 0.0, "init");
  double *b_xi = NULL;
  b_xi = (double *)mem_calloc_1d_double(size*3*3, 0.0, "init");
  double *c_xi = NULL;
  c_xi = (double *)mem_calloc_1d_double(size*3*3, 0.0, "init");
  double *g_xi = NULL;
  g_xi = (double *)mem_calloc_1d_double(size*3, 0.0, "init");

  // et
  double *a_et = NULL;
  a_et = (double *)mem_calloc_1d_double(size*3*3, 0.0, "init");
  double *b_et = NULL;
  b_et = (double *)mem_calloc_1d_double(size*3*3, 0.0, "init");
  double *c_et = NULL;
  c_et = (double *)mem_calloc_1d_double(size*3*3, 0.0, "init");
  double *d_et = NULL;
  d_et = (double *)mem_calloc_1d_double(size*3, 0.0, "init");

  double *xyz = NULL;
  xyz = (double *)mem_calloc_1d_double(size*3, 0.0, "init");

  // solve first layer
  int k=1;
  cal_matrix(gdcurv,k,a_xi,b_xi,c_xi,a_et,b_et,c_et,d_et,volume);
  modify_smooth(gdcurv,k,coef_e_xi,coef_e_et,a_xi,b_xi,c_xi,a_et,b_et,c_et,d_et);
  modify_bdry(a_xi,b_xi,c_xi,a_et,b_et,c_et,d_et,nx,ny,epsilon_x,bdry_x_itype,epsilon_y,bdry_y_itype);
  solve_et_block(a_et,b_et,c_et,d_et,g_xi,nx,ny);
  solve_xi_block(a_xi,b_xi,c_xi,g_xi,xyz,nx,ny);
  assign_coords(gdcurv,xyz,k,epsilon_x,bdry_x_itype,epsilon_y,bdry_y_itype);


  for(int k=1; k<nz; k++)
  {
    cal_smooth_coef(gdcurv,coef,k,coef_e_xi,coef_e_et);
    cal_matrix(gdcurv,k,a_xi,b_xi,c_xi,a_et,b_et,c_et,d_et,volume);
    modify_smooth(gdcurv,k,coef_e_xi,coef_e_et,a_xi,b_xi,c_xi,a_et,b_et,c_et,d_et);
    modify_bdry(a_xi,b_xi,c_xi,a_et,b_et,c_et,d_et,nx,ny,epsilon_x,bdry_x_itype,epsilon_y,bdry_y_itype);
    solve_et_block(a_et,b_et,c_et,d_et,g_xi,nx,ny);
    solve_xi_block(a_xi,b_xi,c_xi,g_xi,xyz,nx,ny);
    assign_coords(gdcurv,xyz,k,epsilon_x,bdry_x_itype,epsilon_y,bdry_y_itype);
    fprintf(stdout,"number of layer is %d\n",k);
    fflush(stdout);
  }

  if(o2i == 1)
  {
    //fprintf(stdout,"hyperbolic method, inner bdry(k=0), outer bdry(nz-1)\n");
    //fprintf(stdout,"we default set read init bdry is inner bdry(k=0)\n");
    //fprintf(stdout,"so if the init bdry is outer bdry actually, must be flip\n");
    flip_coord_z(gdcurv);
  }

  free(coef_e_xi);
  free(coef_e_et);
  free(volume);
  free(a_xi);
  free(b_xi);
  free(c_xi);
  free(g_xi);
  free(a_et);
  free(b_et);
  free(c_et);
  free(d_et);
  free(xyz);

  return 0;
}

int
cal_smooth_coef(gd_t *gdcurv, float coef, int k, float *coef_e_xi, float *coef_e_et)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;

  float S;
  size_t iptr1,iptr2,iptr3,iptr4;
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float xi_len,et_len,zt_len,N_xi,N_et;
  float x_xi_plus, y_xi_plus, z_xi_plus;
  float x_xi_minus,y_xi_minus,z_xi_minus;
  float x_et_plus, y_et_plus, z_et_plus;
  float x_et_minus,y_et_minus,z_et_minus;
  float xi_plus1,xi_minus1,xi_plus2,xi_minus2;
  float et_plus1,et_minus1,et_plus2,et_minus2;
  float d_xi1,d_xi2,d_et1,d_et2;
  float delta_xi,delta_xi_mdfy;
  float delta_et,delta_et_mdfy;
  float len_vec;
  double r_xi[3], r_et[3], vec_n[3];
  float cos_alpha_xi1,cos_alpha_xi2,cos_alpha_xi;
  float cos_alpha_et1,cos_alpha_et2,cos_alpha_et;
  float alpha_xi,alpha_et;

  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;

  S = sqrt((1.0*k)/(nz-1));
  //if(nz<=50) {
  //  S = sqrt((1.0*k)/(nz-1));
  //}

  //if(nz>50) 
  //{
  //  if(k<=50)
  //  {
  //    S = sqrt((1.0*k)/(50));
  //  } else {
  //    S = 1.0 + (k-51)/(nz-51);
  //  }
  //}

  if(k==1)
  {
    int k1=2;
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = (k1-1)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-1)
        iptr2 = (k1-1)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 

        iptr1 = (k1-1)*siz_iz + (j+1)*siz_iy + i; // (i,j+1,k-1)
        iptr2 = (k1-1)*siz_iz + (j-1)*siz_iy + i; // (i,j-1,k-1)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]); 

        iptr1 = (k1-1)*siz_iz + j*siz_iy + i;     // (i,j,k-1)
        iptr2 = (k1-2)*siz_iz + j*siz_iy + i;     // (i,j,k-2)
        x_zt = x3d[iptr1] - x3d[iptr2];
        y_zt = y3d[iptr1] - y3d[iptr2];
        z_zt = z3d[iptr1] - z3d[iptr2];

        xi_len = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        et_len = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        zt_len = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        N_xi = zt_len/xi_len;
        N_et = zt_len/et_len;

        iptr1 = (k1-2)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-2)
        iptr2 = (k1-2)*siz_iz + j*siz_iy + i;     // (i,  j,k-2)
        iptr3 = (k1-2)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-2)
        x_xi_plus = x3d[iptr1] - x3d[iptr2];
        y_xi_plus = y3d[iptr1] - y3d[iptr2];
        z_xi_plus = z3d[iptr1] - z3d[iptr2];
        xi_plus1 = sqrt(pow(x_xi_plus,2) + pow(y_xi_plus,2) + pow(z_xi_plus,2));
        x_xi_minus = x3d[iptr2] - x3d[iptr3];
        y_xi_minus = y3d[iptr2] - y3d[iptr3];
        z_xi_minus = z3d[iptr2] - z3d[iptr3];
        xi_minus1 = sqrt(pow(x_xi_minus,2) + pow(y_xi_minus,2) + pow(z_xi_minus,2));
        iptr1 = (k1-2)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-2)
        iptr2 = (k1-2)*siz_iz + j*siz_iy + i;       // (i,  j,k-2)
        iptr3 = (k1-2)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-2)
        x_et_plus = x3d[iptr1] - x3d[iptr2];
        y_et_plus = y3d[iptr1] - y3d[iptr2];
        z_et_plus = z3d[iptr1] - z3d[iptr2];
        et_plus1 = sqrt(pow(x_et_plus,2) + pow(y_et_plus,2) + pow(z_et_plus,2));
        x_et_minus = x3d[iptr2] - x3d[iptr3];
        y_et_minus = y3d[iptr2] - y3d[iptr3];
        z_et_minus = z3d[iptr2] - z3d[iptr3];
        et_minus1 = sqrt(pow(x_et_minus,2) + pow(y_et_minus,2) + pow(z_et_minus,2));

        iptr1 = (k1-1)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-1)
        iptr2 = (k1-1)*siz_iz + j*siz_iy + i;     // (i,  j,k-1)
        iptr3 = (k1-1)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-1)
        x_xi_plus = x3d[iptr1] - x3d[iptr2];
        y_xi_plus = y3d[iptr1] - y3d[iptr2];
        z_xi_plus = z3d[iptr1] - z3d[iptr2];
        xi_plus2 = sqrt(pow(x_xi_plus,2) + pow(y_xi_plus,2) + pow(z_xi_plus,2));
        x_xi_minus = x3d[iptr2] - x3d[iptr3];
        y_xi_minus = y3d[iptr2] - y3d[iptr3];
        z_xi_minus = z3d[iptr2] - z3d[iptr3];
        xi_minus2 = sqrt(pow(x_xi_minus,2) + pow(y_xi_minus,2) + pow(z_xi_minus,2));
        iptr1 = (k1-1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-1)
        iptr2 = (k1-1)*siz_iz + j*siz_iy + i;       // (i,  j,k-1)
        iptr3 = (k1-1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-1)
        x_et_plus = x3d[iptr1] - x3d[iptr2];
        y_et_plus = y3d[iptr1] - y3d[iptr2];
        z_et_plus = z3d[iptr1] - z3d[iptr2];
        et_plus2 = sqrt(pow(x_et_plus,2) + pow(y_et_plus,2) + pow(z_et_plus,2));
        x_et_minus = x3d[iptr2] - x3d[iptr3];
        y_et_minus = y3d[iptr2] - y3d[iptr3];
        z_et_minus = z3d[iptr2] - z3d[iptr3];
        et_minus2 = sqrt(pow(x_et_minus,2) + pow(y_et_minus,2) + pow(z_et_minus,2));

        d_xi1 = xi_plus1 + xi_minus1;
        d_xi2 = xi_plus2 + xi_minus2;
        d_et1 = et_plus1 + et_minus1;
        d_et2 = et_plus2 + et_minus2;

        delta_xi = d_xi1/d_xi2;
        delta_et = d_et1/d_et2;
        delta_xi_mdfy = fmax(pow(delta_xi,2/S),0.1);
        delta_et_mdfy = fmax(pow(delta_et,2/S),0.1);

        // normalization
        iptr1 = (k1-1)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-1)
        iptr2 = (k1-1)*siz_iz + j*siz_iy + i;     // (i,  j,k-1)
        iptr3 = (k1-1)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-1)
        x_xi_plus = (x3d[iptr1]-x3d[iptr2])/xi_plus2;
        y_xi_plus = (y3d[iptr1]-y3d[iptr2])/xi_plus2;
        z_xi_plus = (z3d[iptr1]-z3d[iptr2])/xi_plus2;
        x_xi_minus = (x3d[iptr3]-x3d[iptr2])/xi_minus2;
        y_xi_minus = (y3d[iptr3]-y3d[iptr2])/xi_minus2;
        z_xi_minus = (z3d[iptr3]-z3d[iptr2])/xi_minus2;
        iptr1 = (k1-1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-1)
        iptr2 = (k1-1)*siz_iz + j*siz_iy + i;       // (i,  j,k-1)
        iptr3 = (k1-1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-1)
        x_et_plus = (x3d[iptr1]-x3d[iptr2])/et_plus2;
        y_et_plus = (y3d[iptr1]-y3d[iptr2])/et_plus2;
        z_et_plus = (z3d[iptr1]-z3d[iptr2])/et_plus2;
        x_et_minus = (x3d[iptr3]-x3d[iptr2])/et_minus2;
        y_et_minus = (y3d[iptr3]-y3d[iptr2])/et_minus2;
        z_et_minus = (z3d[iptr3]-z3d[iptr2])/et_minus2;

        r_xi[0] = x_xi_plus - x_xi_minus; 
        r_xi[1] = y_xi_plus - y_xi_minus; 
        r_xi[2] = z_xi_plus - z_xi_minus; 
        r_et[0] = x_et_plus - x_et_minus; 
        r_et[1] = y_et_plus - y_et_minus; 
        r_et[2] = z_et_plus - z_et_minus; 

        // define vec_n = r_xi X r_et
        cross_product(r_xi,r_et,vec_n);
        len_vec = sqrt(pow(vec_n[0],2) + pow(vec_n[1],2) + pow(vec_n[2],2)); 
        vec_n[0] = vec_n[0]/len_vec;
        vec_n[1] = vec_n[1]/len_vec;
        vec_n[2] = vec_n[2]/len_vec;

        cos_alpha_xi1 = vec_n[0]*x_xi_plus  + vec_n[1]*y_xi_plus  + vec_n[2]*z_xi_plus;
        cos_alpha_xi2 = vec_n[0]*x_xi_minus + vec_n[1]*y_xi_minus + vec_n[2]*z_xi_minus;
        cos_alpha_xi = 0.5*(cos_alpha_xi1 + cos_alpha_xi2);

        cos_alpha_et1 = vec_n[0]*x_et_plus  + vec_n[1]*y_et_plus  + vec_n[2]*z_et_plus;
        cos_alpha_et2 = vec_n[0]*x_et_minus + vec_n[1]*y_et_minus + vec_n[2]*z_et_minus;
        cos_alpha_et = 0.5*(cos_alpha_et1 + cos_alpha_et2);
        if(cos_alpha_xi>=0)
        {
          alpha_xi = 1.0/(1-pow(cos_alpha_xi,2));
        }
        if(cos_alpha_xi<0)
        {
          alpha_xi = 1.0;
        }
        if(cos_alpha_et>=0)
        {
          alpha_et = 1.0/(1-pow(cos_alpha_et,2));
        }
        if(cos_alpha_et<0)
        {
          alpha_et = 1.0;
        }
        if(cos_alpha_xi>1 || cos_alpha_xi<(-1) || cos_alpha_et>1 ||cos_alpha_et<(-1))
        {
          fprintf(stdout,"angle calculation is wrong\n");
          fflush(stdout); exit(1);
        }
        iptr1 = (j-1)*(nx-2) + i-1;
        coef_e_xi[iptr1] = coef*N_xi*S*delta_xi_mdfy*alpha_xi;
        coef_e_et[iptr1] = coef*N_et*S*delta_et_mdfy*alpha_et;
      }
    }
  }


  if(k>1)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = (k-1)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 

        iptr1 = (k-1)*siz_iz + (j+1)*siz_iy + i; // (i,j+1,k-1)
        iptr2 = (k-1)*siz_iz + (j-1)*siz_iy + i; // (i,j-1,k-1)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]); 

        iptr1 = (k-1)*siz_iz + j*siz_iy + i;     // (i,j,k-1)
        iptr2 = (k-2)*siz_iz + j*siz_iy + i;     // (i,j,k-2)
        x_zt = x3d[iptr1] - x3d[iptr2];
        y_zt = y3d[iptr1] - y3d[iptr2];
        z_zt = z3d[iptr1] - z3d[iptr2];

        xi_len = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        et_len = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        zt_len = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        N_xi = zt_len/xi_len;
        N_et = zt_len/et_len;

        iptr1 = (k-2)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-2)
        iptr2 = (k-2)*siz_iz + j*siz_iy + i;     // (i,  j,k-2)
        iptr3 = (k-2)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-2)
        x_xi_plus = x3d[iptr1] - x3d[iptr2];
        y_xi_plus = y3d[iptr1] - y3d[iptr2];
        z_xi_plus = z3d[iptr1] - z3d[iptr2];
        xi_plus1 = sqrt(pow(x_xi_plus,2) + pow(y_xi_plus,2) + pow(z_xi_plus,2));
        x_xi_minus = x3d[iptr2] - x3d[iptr3];
        y_xi_minus = y3d[iptr2] - y3d[iptr3];
        z_xi_minus = z3d[iptr2] - z3d[iptr3];
        xi_minus1 = sqrt(pow(x_xi_minus,2) + pow(y_xi_minus,2) + pow(z_xi_minus,2));
        iptr1 = (k-2)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-2)
        iptr2 = (k-2)*siz_iz + j*siz_iy + i;       // (i,  j,k-2)
        iptr3 = (k-2)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-2)
        x_et_plus = x3d[iptr1] - x3d[iptr2];
        y_et_plus = y3d[iptr1] - y3d[iptr2];
        z_et_plus = z3d[iptr1] - z3d[iptr2];
        et_plus1 = sqrt(pow(x_et_plus,2) + pow(y_et_plus,2) + pow(z_et_plus,2));
        x_et_minus = x3d[iptr2] - x3d[iptr3];
        y_et_minus = y3d[iptr2] - y3d[iptr3];
        z_et_minus = z3d[iptr2] - z3d[iptr3];
        et_minus1 = sqrt(pow(x_et_minus,2) + pow(y_et_minus,2) + pow(z_et_minus,2));

        iptr1 = (k-1)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + i;     // (i,  j,k-1)
        iptr3 = (k-1)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-1)
        x_xi_plus = x3d[iptr1] - x3d[iptr2];
        y_xi_plus = y3d[iptr1] - y3d[iptr2];
        z_xi_plus = z3d[iptr1] - z3d[iptr2];
        xi_plus2 = sqrt(pow(x_xi_plus,2) + pow(y_xi_plus,2) + pow(z_xi_plus,2));
        x_xi_minus = x3d[iptr2] - x3d[iptr3];
        y_xi_minus = y3d[iptr2] - y3d[iptr3];
        z_xi_minus = z3d[iptr2] - z3d[iptr3];
        xi_minus2 = sqrt(pow(x_xi_minus,2) + pow(y_xi_minus,2) + pow(z_xi_minus,2));
        iptr1 = (k-1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + i;       // (i,  j,k-1)
        iptr3 = (k-1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-1)
        x_et_plus = x3d[iptr1] - x3d[iptr2];
        y_et_plus = y3d[iptr1] - y3d[iptr2];
        z_et_plus = z3d[iptr1] - z3d[iptr2];
        et_plus2 = sqrt(pow(x_et_plus,2) + pow(y_et_plus,2) + pow(z_et_plus,2));
        x_et_minus = x3d[iptr2] - x3d[iptr3];
        y_et_minus = y3d[iptr2] - y3d[iptr3];
        z_et_minus = z3d[iptr2] - z3d[iptr3];
        et_minus2 = sqrt(pow(x_et_minus,2) + pow(y_et_minus,2) + pow(z_et_minus,2));

        d_xi1 = xi_plus1 + xi_minus1;
        d_xi2 = xi_plus2 + xi_minus2;
        d_et1 = et_plus1 + et_minus1;
        d_et2 = et_plus2 + et_minus2;

        delta_xi = d_xi1/d_xi2;
        delta_et = d_et1/d_et2;
        delta_xi_mdfy = fmax(pow(delta_xi,2/S),0.1);
        delta_et_mdfy = fmax(pow(delta_et,2/S),0.1);

        // normalization
        iptr1 = (k-1)*siz_iz + j*siz_iy + i+1;   // (i+1,j,k-1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + i;     // (i,  j,k-1)
        iptr3 = (k-1)*siz_iz + j*siz_iy + i-1;   // (i-1,j,k-1)
        x_xi_plus = (x3d[iptr1]-x3d[iptr2])/xi_plus2;
        y_xi_plus = (y3d[iptr1]-y3d[iptr2])/xi_plus2;
        z_xi_plus = (z3d[iptr1]-z3d[iptr2])/xi_plus2;
        x_xi_minus = (x3d[iptr3]-x3d[iptr2])/xi_minus2;
        y_xi_minus = (y3d[iptr3]-y3d[iptr2])/xi_minus2;
        z_xi_minus = (z3d[iptr3]-z3d[iptr2])/xi_minus2;
        iptr1 = (k-1)*siz_iz + (j+1)*siz_iy + i;   // (i,j+1,k-1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + i;       // (i,  j,k-1)
        iptr3 = (k-1)*siz_iz + (j-1)*siz_iy + i;   // (i,j-1,k-1)
        x_et_plus = (x3d[iptr1]-x3d[iptr2])/et_plus2;
        y_et_plus = (y3d[iptr1]-y3d[iptr2])/et_plus2;
        z_et_plus = (z3d[iptr1]-z3d[iptr2])/et_plus2;
        x_et_minus = (x3d[iptr3]-x3d[iptr2])/et_minus2;
        y_et_minus = (y3d[iptr3]-y3d[iptr2])/et_minus2;
        z_et_minus = (z3d[iptr3]-z3d[iptr2])/et_minus2;

        r_xi[0] = x_xi_plus - x_xi_minus; 
        r_xi[1] = y_xi_plus - y_xi_minus; 
        r_xi[2] = z_xi_plus - z_xi_minus; 
        r_et[0] = x_et_plus - x_et_minus; 
        r_et[1] = y_et_plus - y_et_minus; 
        r_et[2] = z_et_plus - z_et_minus; 

        // define vec_n = r_xi X r_et
        cross_product(r_xi,r_et,vec_n);
        len_vec = sqrt(pow(vec_n[0],2) + pow(vec_n[1],2) + pow(vec_n[2],2)); 
        vec_n[0] = vec_n[0]/len_vec;
        vec_n[1] = vec_n[1]/len_vec;
        vec_n[2] = vec_n[2]/len_vec;

        cos_alpha_xi1 = vec_n[0]*x_xi_plus  + vec_n[1]*y_xi_plus  + vec_n[2]*z_xi_plus;
        cos_alpha_xi2 = vec_n[0]*x_xi_minus + vec_n[1]*y_xi_minus + vec_n[2]*z_xi_minus;
        cos_alpha_xi = 0.5*(cos_alpha_xi1 + cos_alpha_xi2);

        cos_alpha_et1 = vec_n[0]*x_et_plus  + vec_n[1]*y_et_plus  + vec_n[2]*z_et_plus;
        cos_alpha_et2 = vec_n[0]*x_et_minus + vec_n[1]*y_et_minus + vec_n[2]*z_et_minus;
        cos_alpha_et = 0.5*(cos_alpha_et1 + cos_alpha_et2);
        if(cos_alpha_xi>=0)
        {
          alpha_xi = 1.0/(1-pow(cos_alpha_xi,2));
        }
        if(cos_alpha_xi<0)
        {
          alpha_xi = 1.0;
        }
        if(cos_alpha_et>=0)
        {
          alpha_et = 1.0/(1-pow(cos_alpha_et,2));
        }
        if(cos_alpha_et<0)
        {
          alpha_et = 1.0;
        }
        if(cos_alpha_xi>1 || cos_alpha_xi<(-1) || cos_alpha_et>1 ||cos_alpha_et<(-1))
        {
          fprintf(stdout,"angle calculation is wrong\n");
          fflush(stdout); exit(1);
        }
        iptr1 = (j-1)*(nx-2) + i-1;
        coef_e_xi[iptr1] = coef*N_xi*S*delta_xi_mdfy*alpha_xi;
        coef_e_et[iptr1] = coef*N_et*S*delta_et_mdfy*alpha_et;
      }
    }
  }

  return 0;
}

int 
cal_matrix(gd_t *gdcurv, int k, double *a_xi, double *b_xi, double *c_xi, 
           double *a_et, double *b_et, double *c_et, double *d_et, float *volume)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *step = gdcurv->step;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  
  double A[3][3],B[3][3],C[3][3];
  double mat1[3][3], mat2[3][3], mat3[3][3];
  double vec_d[3], vec_g[3];
  size_t iptr,iptr1,iptr2;
  double r_xi0[3],r_et0[3],r_zt0[3];
  double sigma[3],omega[3],tau[3];
  double area;


  if(k==1)
  {
    for(int j=1; j<ny-1; j++) 
    {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = (k-1)*siz_iz + j*siz_iy + i+1;
        iptr2 = (k-1)*siz_iz + j*siz_iy + i-1;
        r_xi0[0] = 0.5*(x3d[iptr1] - x3d[iptr2]);
        r_xi0[1] = 0.5*(y3d[iptr1] - y3d[iptr2]);
        r_xi0[2] = 0.5*(z3d[iptr1] - z3d[iptr2]);

        iptr1 = (k-1)*siz_iz + (j+1)*siz_iy + i;
        iptr2 = (k-1)*siz_iz + (j-1)*siz_iy + i;
        r_et0[0] = 0.5*(x3d[iptr1] - x3d[iptr2]);
        r_et0[1] = 0.5*(y3d[iptr1] - y3d[iptr2]);
        r_et0[2] = 0.5*(z3d[iptr1] - z3d[iptr2]);

        // define sigma = r_xi0 X r_et0
        cross_product(r_xi0,r_et0,sigma);

        C[0][0] = r_xi0[0]+1e-7; C[0][1] = r_xi0[1];      C[0][2] = r_xi0[2];
        C[1][0] = r_et0[0]     ; C[1][1] = r_et0[1]+1e-7; C[1][2] = r_et0[2];
        C[2][0] = sigma[0]     ; C[2][1] = sigma[1];      C[2][2] = sigma[2]+1e-7;

        area = sqrt(pow(sigma[0],2) + pow(sigma[1],2) + pow(sigma[2],2));
        // area->volume
        // volume(j,i,1) = V1 volume(j,i,2) = V0
        iptr = j*nx + i; 
        volume[iptr] = area * step[k-1];
        volume[iptr+siz_iz] = volume[iptr];

        vec_g[0] = 0; vec_g[1] = 0; vec_g[2] = volume[iptr+siz_iz];

        // r_zt0 = inv(C) * g
        mat_invert3x3(C);
        mat_mul3x1(C,vec_g,r_zt0);

        // define omega = r_zt0 X r_xi0
        cross_product(r_zt0,r_xi0,omega);

        // define tau = r_et0 X r_zt0
        cross_product(r_et0,r_zt0,tau);

        A[0][0] = r_zt0[0]; A[0][1] = r_zt0[1]; A[0][2] = r_zt0[2];
        A[1][0] = 0       ; A[1][1] = 0       ; A[1][2] = 0       ;
        A[2][0] = tau[0]  ; A[2][1] = tau[1]  ; A[2][2] = tau[2];

        B[0][0] = 0;        B[0][1] = 0       ; B[0][2] = 0;
        B[1][0] = r_zt0[0]; B[1][1] = r_zt0[1]; B[1][2] = r_zt0[2];
        B[2][0] = omega[0]; B[2][1] = omega[1]; B[2][2] = omega[2];

        mat_mul3x3(C,A,mat1); //mat1 = inv(C) * A
        mat_mul3x3(C,B,mat2); //mat2 = inv(C) * B
        mat_iden3x3(mat3);
        // NOTE: now volume is volume(j,i,1) = V1
        vec_g[0] = 0; vec_g[1] = 0; vec_g[2] = volume[iptr];
        mat_mul3x1(C,vec_g,vec_d);
        iptr1 = ((j-1)*(nx-2) + (i-1))*3*3;
        iptr2 = ((j-1)*(nx-2) + (i-1))*3;
        for(int ii=0; ii<3; ii++) {
          for(int jj=0; jj<3; jj++) {
            a_xi[iptr1+3*ii+jj] = -0.5*mat1[ii][jj];
            b_xi[iptr1+3*ii+jj] = mat3[ii][jj];
            c_xi[iptr1+3*ii+jj] = 0.5*mat1[ii][jj];

            a_et[iptr1+3*ii+jj] = -0.5*mat2[ii][jj];
            b_et[iptr1+3*ii+jj] = mat3[ii][jj];
            c_et[iptr1+3*ii+jj] = 0.5*mat2[ii][jj];
          }
          d_et[iptr2+ii] = vec_d[ii];
        }
      }
    }
  } else {
    for(int j=1; j<ny-1; j++) 
    {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = (k-1)*siz_iz + j*siz_iy + i+1;
        iptr2 = (k-1)*siz_iz + j*siz_iy + i-1;
        r_xi0[0] = 0.5*(x3d[iptr1] - x3d[iptr2]);
        r_xi0[1] = 0.5*(y3d[iptr1] - y3d[iptr2]);
        r_xi0[2] = 0.5*(z3d[iptr1] - z3d[iptr2]);

        iptr1 = (k-1)*siz_iz + (j+1)*siz_iy + i;
        iptr2 = (k-1)*siz_iz + (j-1)*siz_iy + i;
        r_et0[0] = 0.5*(x3d[iptr1] - x3d[iptr2]);
        r_et0[1] = 0.5*(y3d[iptr1] - y3d[iptr2]);
        r_et0[2] = 0.5*(z3d[iptr1] - z3d[iptr2]);

        // define sigma = r_xi0 X r_et0
        cross_product(r_xi0,r_et0,sigma);

        // add damping factor, maybe inv(C) is singular
        C[0][0] = r_xi0[0]+1e-7; C[0][1] = r_xi0[1];      C[0][2] = r_xi0[2];
        C[1][0] = r_et0[0]     ; C[1][1] = r_et0[1]+1e-7; C[1][2] = r_et0[2];
        C[2][0] = sigma[0]     ; C[2][1] = sigma[1];      C[2][2] = sigma[2]+1e-7;

        area = sqrt(pow(sigma[0],2) + pow(sigma[1],2) + pow(sigma[2],2));
        // area->volume
        // volume(j,i,1) = V1 volume(j,i,2) = V0
        iptr = j*nx + i; 
        volume[iptr+siz_iz] = volume[iptr];
        volume[iptr] = area * step[k-1];

        vec_g[0] = 0; vec_g[1] = 0; vec_g[2] = volume[iptr+siz_iz];

        // r_zt0 = inv(C) * g
        mat_invert3x3(C);
        mat_mul3x1(C,vec_g,r_zt0);

        // define omega = r_zt0 X r_xi0
        cross_product(r_zt0,r_xi0,omega);

        // define tau = r_et0 X r_zt0
        cross_product(r_et0,r_zt0,tau);

        A[0][0] = r_zt0[0]; A[0][1] = r_zt0[1]; A[0][2] = r_zt0[2];
        A[1][0] = 0       ; A[1][1] = 0       ; A[1][2] = 0       ;
        A[2][0] = tau[0]  ; A[2][1] = tau[1]  ; A[2][2] = tau[2];

        B[0][0] = 0;        B[0][1] = 0       ; B[0][2] = 0;
        B[1][0] = r_zt0[0]; B[1][1] = r_zt0[1]; B[1][2] = r_zt0[2];
        B[2][0] = omega[0]; B[2][1] = omega[1]; B[2][2] = omega[2];

        mat_mul3x3(C,A,mat1); //mat1 = inv(C) * A
        mat_mul3x3(C,B,mat2); //mat2 = inv(C) * B
        mat_iden3x3(mat3);
        // NOTE: now volume is volume(j,i,1) = V1
        vec_g[0] = 0; vec_g[1] = 0; vec_g[2] = volume[iptr];
        mat_mul3x1(C,vec_g,vec_d);
        iptr1 = ((j-1)*(nx-2) + (i-1))*3*3;
        iptr2 = ((j-1)*(nx-2) + (i-1))*3;
        for(int ii=0; ii<3; ii++) {
          for(int jj=0; jj<3; jj++) {
            a_xi[iptr1+3*ii+jj] = -0.5*mat1[ii][jj];
            b_xi[iptr1+3*ii+jj] = mat3[ii][jj];
            c_xi[iptr1+3*ii+jj] = 0.5*mat1[ii][jj];

            a_et[iptr1+3*ii+jj] = -0.5*mat2[ii][jj];
            b_et[iptr1+3*ii+jj] = mat3[ii][jj];
            c_et[iptr1+3*ii+jj] = 0.5*mat2[ii][jj];
          }
          d_et[iptr2+ii] = vec_d[ii];
        }
      }
    }
  }

  return 0;
}

int
modify_smooth(gd_t *gdcurv, int k, float *coef_e_xi, float *coef_e_et, 
              double *a_xi, double *b_xi, double *c_xi, double *a_et, 
              double *b_et, double *c_et, double *d_et)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;

  double vec1[3],vec2[3],vec3[3],vec4[3],vec5[3];
  float coef_i_xi,coef_i_et;
  size_t iptr1, iptr2, iptr3, iptr4, iptr5;
  double mat[3][3];
  mat_iden3x3(mat);
  for(int j=1; j<ny-1; j++)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*siz_iz + j*siz_iy + i;    //(i,j,k-1)
      iptr2 = (k-1)*siz_iz + j*siz_iy + i+1;  //(i+1,j,k-1)
      iptr3 = (k-1)*siz_iz + j*siz_iy + i-1;  //(i-1,j,k-1)
      iptr4 = (k-1)*siz_iz + (j+1)*siz_iy + i;  //(i,j+1,k-1)
      iptr5 = (k-1)*siz_iz + (j-1)*siz_iy + i;  //(i,j-1,k-1)
      vec1[0] = x3d[iptr1]; vec1[1] = y3d[iptr1]; vec1[2] = z3d[iptr1];
      vec2[0] = x3d[iptr2]; vec2[1] = y3d[iptr2]; vec2[2] = z3d[iptr2];
      vec3[0] = x3d[iptr3]; vec3[1] = y3d[iptr3]; vec3[2] = z3d[iptr3];
      vec4[0] = x3d[iptr4]; vec4[1] = y3d[iptr4]; vec4[2] = z3d[iptr4];
      vec5[0] = x3d[iptr5]; vec5[1] = y3d[iptr5]; vec5[2] = z3d[iptr5];

      iptr1 = ((j-1)*(nx-2)+(i-1))*3*3;
      iptr2 = ((j-1)*(nx-2)+(i-1))*3;
      iptr3 = (j-1)*(nx-2)+(i-1);
      
      coef_i_xi = 2*coef_e_xi[iptr3];
      coef_i_et = 2*coef_e_et[iptr3];

      for(int ii=0; ii<3; ii++) {
        for(int jj=0; jj<3; jj++) {
          a_xi[iptr1+3*ii+jj] = a_xi[iptr1+3*ii+jj] - coef_i_xi*mat[ii][jj];
          b_xi[iptr1+3*ii+jj] = b_xi[iptr1+3*ii+jj] + 2*coef_i_xi*mat[ii][jj];
          c_xi[iptr1+3*ii+jj] = c_xi[iptr1+3*ii+jj] - coef_i_xi*mat[ii][jj];

          a_et[iptr1+3*ii+jj] = a_et[iptr1+3*ii+jj] - coef_i_et*mat[ii][jj];
          b_et[iptr1+3*ii+jj] = b_et[iptr1+3*ii+jj] + 2*coef_i_et*mat[ii][jj];
          c_et[iptr1+3*ii+jj] = c_et[iptr1+3*ii+jj] - coef_i_et*mat[ii][jj];
        }
        d_et[iptr2+ii] = d_et[iptr2+ii] + coef_e_xi[iptr3]*(vec2[ii]+vec3[ii]-2*vec1[ii])
                       + coef_e_et[iptr3]*(vec4[ii]+vec5[ii]-2*vec1[ii]);
      }
    }
  }

  return 0;
}

int
modify_bdry(double *a_xi, double *b_xi, double *c_xi, double *a_et, 
            double *b_et, double *c_et, double *d_et, int nx, int ny,
            float epsilon_x, int bdry_x_itype,
            float epsilon_y, int bdry_y_itype)
{
  size_t iptr1,iptr2;
  // float boundry
  if(bdry_x_itype == 1)
  {
    for(int j=1; j<ny-1; j++)
    {
      iptr1 = ((j-1)*(nx-2)+0)*3*3;
      iptr2 = ((j-1)*(nx-2)+(nx-3))*3*3;
      for(int ii=0; ii<3; ii++) {
        for(int jj=0; jj<3; jj++) {
          // modify i=0
          b_xi[iptr1+ii*3+jj] = b_xi[iptr1+ii*3+jj] + (1+epsilon_x)*a_xi[iptr1+ii*3+jj];
          c_xi[iptr1+ii*3+jj] = c_xi[iptr1+ii*3+jj] - epsilon_x*a_xi[iptr1+ii*3+jj];
          // modify i=nx-3
          b_xi[iptr2+ii*3+jj] = b_xi[iptr2+ii*3+jj] + (1+epsilon_x)*c_xi[iptr2+ii*3+jj];
          a_xi[iptr2+ii*3+jj] = a_xi[iptr2+ii*3+jj] - epsilon_x*c_xi[iptr2+ii*3+jj];
        }
      }
    }
  }
  // cartesian boundary dx=0
  if(bdry_x_itype == 2)
  {
    for(int j=1; j<ny-1; j++)
    {
      iptr1 = ((j-1)*(nx-2)+0)*3*3;
      iptr2 = ((j-1)*(nx-2)+(nx-3))*3*3;
      for(int ii=0; ii<3; ii++) {
        for(int jj=0; jj<3; jj++) {
          // only modify second and third column 
          // due to dx=0, so first colume equal to 0
          if(jj == 0) continue;
          // modify i=0
          b_xi[iptr1+ii*3+jj] = b_xi[iptr1+ii*3+jj] + a_xi[iptr1+ii*3+jj];
          // modify i=nx-3
          b_xi[iptr2+ii*3+jj] = b_xi[iptr2+ii*3+jj] + c_xi[iptr2+ii*3+jj];
        }
      }
    }
  }

  if(bdry_y_itype == 1)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (0*(nx-2)+(i-1))*3*3;
      iptr2 = ((ny-3)*(nx-2)+(i-1))*3*3;
      for(int ii=0; ii<3; ii++) {
        for(int jj=0; jj<3; jj++) {
          // modify j=0
          b_et[iptr1+ii*3+jj] = b_et[iptr1+ii*3+jj] + (1+epsilon_y)*a_et[iptr1+ii*3+jj];
          c_et[iptr1+ii*3+jj] = c_et[iptr1+ii*3+jj] - epsilon_y*a_et[iptr1+ii*3+jj];
          // modify j=ny-3
          b_et[iptr2+ii*3+jj] = b_et[iptr2+ii*3+jj] + (1+epsilon_y)*c_et[iptr2+ii*3+jj];
          a_et[iptr2+ii*3+jj] = a_et[iptr2+ii*3+jj] - epsilon_y*c_et[iptr2+ii*3+jj];
        }
      }
    }
  }

  // cartesian boundary dy=0
  if(bdry_y_itype == 2)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (0*(nx-2)+i-1)*3*3;
      iptr2 = ((ny-3)*(nx-2)+i-1)*3*3;
      for(int ii=0; ii<3; ii++) {
        for(int jj=0; jj<3; jj++) {
          // only modify first and third column 
          // due to dy=0, so second colume equal to 0
          if(jj == 1) continue;
          // modify j=0
          b_et[iptr1+ii*3+jj] = b_et[iptr1+ii*3+jj] + a_et[iptr1+ii*3+jj];
          // modify j=ny-3
          b_et[iptr2+ii*3+jj] = b_et[iptr2+ii*3+jj] + c_et[iptr2+ii*3+jj];
        }
      }
    }
  }

  return 0;
}

int
assign_coords(gd_t *gdcurv, double *xyz, int k, float epsilon_x, 
              int bdry_x_itype, float epsilon_y, int bdry_y_itype)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;

  size_t iptr,iptr1,iptr2;
  size_t iptr3,iptr4,iptr5;
  for(int j=1; j<ny-1; j++)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr  =  k*siz_iz + j*siz_iy + i;
      iptr1 = (k-1)*siz_iz + j*siz_iy  + i;
      iptr2 = ((j-1)*(nx-2) + (i-1))*3;
      x3d[iptr] = x3d[iptr1] + xyz[iptr2];
      y3d[iptr] = y3d[iptr1] + xyz[iptr2+1];
      z3d[iptr] = z3d[iptr1] + xyz[iptr2+2];
    }
  }

  // floating boundary
  if(bdry_x_itype == 1)
  {
    for(int j=0; j<ny; j++)
    {
      iptr  = k*siz_iz + j*siz_iy + 0;       // (0,j,k)
      iptr1 = (k-1)*siz_iz + j*siz_iy + 0;   // (0,j,k-1)
      iptr2 = k*siz_iz + j*siz_iy +1;        // (1,j,k)
      iptr3 = (k-1)*siz_iz + j*siz_iy + 1;   // (1,j,k-1)
      iptr4 = k*siz_iz + j*siz_iy +2;        // (2,j,k)
      iptr5 = (k-1)*siz_iz + j*siz_iy + 2;   // (2,j,k-1)

      x3d[iptr] = x3d[iptr1] + (1+epsilon_x)*(x3d[iptr2]-x3d[iptr3])
                 -epsilon_x*(x3d[iptr4]-x3d[iptr5]);
      y3d[iptr] = y3d[iptr1] + (1+epsilon_x)*(y3d[iptr2]-y3d[iptr3])
                 -epsilon_x*(y3d[iptr4]-y3d[iptr5]);
      z3d[iptr] = z3d[iptr1] + (1+epsilon_x)*(z3d[iptr2]-z3d[iptr3])
                 -epsilon_x*(z3d[iptr4]-z3d[iptr5]);

      iptr  = k*siz_iz + j*siz_iy + nx-1;       // (nx-1,j,k)
      iptr1 = (k-1)*siz_iz + j*siz_iy + nx-1;   // (nx-1,j,k-1)
      iptr2 = k*siz_iz + j*siz_iy + nx-2;       // (nx-2,j,k)
      iptr3 = (k-1)*siz_iz + j*siz_iy + nx-2;   // (nx-2,j,k-1)
      iptr4 = k*siz_iz + j*siz_iy + nx-3;       // (nx-3,j,k)
      iptr5 = (k-1)*siz_iz + j*siz_iy + nx-3;   // (nx-3,j,k-1)

      x3d[iptr] = x3d[iptr1] + (1+epsilon_x)*(x3d[iptr2]-x3d[iptr3])
                 -epsilon_x*(x3d[iptr4]-x3d[iptr5]);
      y3d[iptr] = y3d[iptr1] + (1+epsilon_x)*(y3d[iptr2]-y3d[iptr3])
                 -epsilon_x*(y3d[iptr4]-y3d[iptr5]);
      z3d[iptr] = z3d[iptr1] + (1+epsilon_x)*(z3d[iptr2]-z3d[iptr3])
                 -epsilon_x*(z3d[iptr4]-z3d[iptr5]);
    }
  }
  // cartesian boundary
  if(bdry_x_itype == 2)
  {
    for(int j=0; j<ny; j++)
    {
      iptr  = k*siz_iz + j*siz_iy + 0;       // (0,j,k)
      iptr1 = (k-1)*siz_iz + j*siz_iy + 0;   // (0,j,k-1)
      iptr2 = k*siz_iz + j*siz_iy +1;        // (1,j,k)
      iptr3 = (k-1)*siz_iz + j*siz_iy + 1;   // (1,j,k-1)

      x3d[iptr] = x3d[iptr1];
      y3d[iptr] = y3d[iptr1] + (y3d[iptr2]-y3d[iptr3]);
      z3d[iptr] = z3d[iptr1] + (z3d[iptr2]-z3d[iptr3]);

      iptr  = k*siz_iz + j*siz_iy + nx-1;       // (nx-1,j,k)
      iptr1 = (k-1)*siz_iz + j*siz_iy + nx-1;   // (nx-1,j,k-1)
      iptr2 = k*siz_iz + j*siz_iy + nx-2;       // (nx-2,j,k)
      iptr3 = (k-1)*siz_iz + j*siz_iy + nx-2;   // (nx-2,j,k-1)

      x3d[iptr] = x3d[iptr1];
      y3d[iptr] = y3d[iptr1] + (y3d[iptr2]-y3d[iptr3]);
      z3d[iptr] = z3d[iptr1] + (z3d[iptr2]-z3d[iptr3]);
    }
  }

  if(bdry_y_itype == 1)
  {
    for(int i=0; i<nx; i++)
    {
      iptr  = k*siz_iz + 0*siz_iy + i;       // (i,0,k)
      iptr1 = (k-1)*siz_iz + 0*siz_iy + i;   // (i,0,k-1)
      iptr2 = k*siz_iz + 1*siz_iy +i;        // (i,1,k)
      iptr3 = (k-1)*siz_iz + 1*siz_iy + i;   // (i,1,k-1)
      iptr4 = k*siz_iz + 2*siz_iy + i;       // (i,2,k)
      iptr5 = (k-1)*siz_iz + 2*siz_iy + i;   // (i,2,k-1)

      x3d[iptr] = x3d[iptr1] + (1+epsilon_y)*(x3d[iptr2]-x3d[iptr3])
                 -epsilon_y*(x3d[iptr4]-x3d[iptr5]);
      y3d[iptr] = y3d[iptr1] + (1+epsilon_y)*(y3d[iptr2]-y3d[iptr3])
                 -epsilon_y*(y3d[iptr4]-y3d[iptr5]);
      z3d[iptr] = z3d[iptr1] + (1+epsilon_y)*(z3d[iptr2]-z3d[iptr3])
                 -epsilon_y*(z3d[iptr4]-z3d[iptr5]);

      iptr  = k*siz_iz + (ny-1)*siz_iy + i;       // (i,ny-1,k)
      iptr1 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   // (i,ny-1,k-1)
      iptr2 = k*siz_iz + (ny-2)*siz_iy + i;       // (i,ny-2,k)
      iptr3 = (k-1)*siz_iz + (ny-2)*siz_iy + i;   // (i,ny-2,k-1)
      iptr4 = k*siz_iz + (ny-3)*siz_iy + i;       // (i,ny-3,k)
      iptr5 = (k-1)*siz_iz + (ny-3)*siz_iy + i;   // (i,ny-3,k-1)

      x3d[iptr] = x3d[iptr1] + (1+epsilon_y)*(x3d[iptr2]-x3d[iptr3])
                 -epsilon_y*(x3d[iptr4]-x3d[iptr5]);
      y3d[iptr] = y3d[iptr1] + (1+epsilon_y)*(y3d[iptr2]-y3d[iptr3])
                 -epsilon_y*(y3d[iptr4]-y3d[iptr5]);
      z3d[iptr] = z3d[iptr1] + (1+epsilon_y)*(z3d[iptr2]-z3d[iptr3])
                 -epsilon_y*(z3d[iptr4]-z3d[iptr5]);
    }
  }

  if(bdry_y_itype == 2)
  {
    for(int i=0; i<nx; i++)
    {
      iptr  = k*siz_iz + 0*siz_iy + i;       // (i,0,k)
      iptr1 = (k-1)*siz_iz + 0*siz_iy + i;   // (i,0,k-1)
      iptr2 = k*siz_iz + 1*siz_iy +i;        // (i,1,k)
      iptr3 = (k-1)*siz_iz + 1*siz_iy + i;   // (i,1,k-1)

      x3d[iptr] = x3d[iptr1] + (x3d[iptr2]-x3d[iptr3]);
      y3d[iptr] = y3d[iptr1];
      z3d[iptr] = z3d[iptr1] + (z3d[iptr2]-z3d[iptr3]);

      iptr  = k*siz_iz + (ny-1)*siz_iy + i;       // (i,ny-1,k)
      iptr1 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   // (i,ny-1,k-1)
      iptr2 = k*siz_iz + (ny-2)*siz_iy + i;       // (i,ny-2,k)
      iptr3 = (k-1)*siz_iz + (ny-2)*siz_iy + i;   // (i,ny-2,k-1)

      x3d[iptr] = x3d[iptr1] + (x3d[iptr2]-x3d[iptr3]);
      y3d[iptr] = y3d[iptr1];
      z3d[iptr] = z3d[iptr1] + (z3d[iptr2]-z3d[iptr3]);
    }
  }

  return 0;
}

int
solve_et_block(double *a_et, double *b_et, double *c_et,
               double *d_et, double *g_xi, int nx, int ny)
{
  size_t iptr1, iptr2;
  int n = ny-2; 
  double *a = NULL;
  double *b = NULL;
  double *c = NULL;
  double *d = NULL;
  double *g = NULL;
  a = (double *)mem_calloc_1d_double(n*3*3, 0.0, "init");
  b = (double *)mem_calloc_1d_double(n*3*3, 0.0, "init");
  c = (double *)mem_calloc_1d_double(n*3*3, 0.0, "init");
  d = (double *)mem_calloc_1d_double(n*3, 0.0, "init");
  g = (double *)mem_calloc_1d_double(n*3, 0.0, "init");
  // temp var
  double *D = NULL;
  double *y = NULL;
  D = (double *)mem_calloc_1d_double(n*3*3, 0.0, "init");
  y = (double *)mem_calloc_1d_double(n*3, 0.0, "init");
  for(int i=0; i<nx-2; i++)
  {
    // get i=const matrix
    for(int j=0; j<ny-2; j++)
    {
      iptr1 = (j*(nx-2) + i)*3*3;
      iptr2 = (j*(nx-2) + i)*3;
      for(int ii=0; ii<3; ii++)
      {
        for(int jj=0; jj<3; jj++)
        {
          a[j*3*3+3*ii+jj] = a_et[iptr1+3*ii+jj];
          b[j*3*3+3*ii+jj] = b_et[iptr1+3*ii+jj];
          c[j*3*3+3*ii+jj] = c_et[iptr1+3*ii+jj];
        }
        d[j*3+ii] = d_et[iptr2+ii];
      }
    }

    thomas_block(n,a,b,c,d,g,D,y);

    for(int j=0; j<ny-2; j++)
    {
      iptr2 = (j*(nx-2) + i)*3;
      for(int ii=0; ii<3; ii++)
      {
        g_xi[iptr2+ii] = g[j*3+ii];
      }
    }

  }

  free(a);
  free(b);
  free(c);
  free(d);
  free(g);
  free(D);
  free(y);

  return 0;
}

int
solve_xi_block(double *a_xi, double *b_xi, double *c_xi,
               double *g_xi, double *xyz, int nx, int ny)
{
  size_t iptr1, iptr2;
  int n = nx-2; 
  double *a;
  double *b;
  double *c;
  double *g;
  double *cord = NULL;
  cord = (double *)mem_calloc_1d_double(n*3, 0.0, "init");
  // temp var
  double *D = NULL;
  double *y = NULL;
  D = (double *)mem_calloc_1d_double(n*3*3, 0.0, "init");
  y = (double *)mem_calloc_1d_double(n*3, 0.0, "init");
  for(int j=0; j<ny-2; j++)
  {
    // get j=const matrix
    iptr1 = (j*(nx-2)+0)*3*3;
    iptr2 = (j*(nx-2)+0)*3;
    a = a_xi + iptr1; 
    b = b_xi + iptr1; 
    c = c_xi + iptr1; 
    g = g_xi + iptr2;

    thomas_block(n,a,b,c,g,cord,D,y);

    for(int i=0; i<nx-2; i++)
    {
      iptr2 = (j*(nx-2) + i)*3;
      for(int ii=0; ii<3; ii++)
      {
        xyz[iptr2+ii] = cord[i*3+ii];
      }
    }

  }

  free(cord);
  free(D);
  free(y);

  return 0;
}
