#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lib_mem.h"
#include "gd_t.h"
#include "constants.h"
#include "algebra.h"
#include "io_funcs.h"

int
init_gdcurv(gd_t *gdcurv, int nx, int nz)
{
  gdcurv->nx = nx;
  gdcurv->nz = nz;
  //2 dimension, x and z
  gdcurv->ncmp = CONST_NDIM; 
  gdcurv->siz_iz   = gdcurv->nx;
  gdcurv->siz_icmp = gdcurv->nx * gdcurv->nz;
  // malloc grid space 
  gdcurv->v3d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }

  // position of each v3d
  size_t *cmp_pos = (size_t *) mem_calloc_1d_sizet(gdcurv->ncmp,
                                                         0,
                                                         "gd_curv_init");
  
  // name of each v3d
  char **cmp_name = (char **) mem_malloc_2l_char(gdcurv->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "gd_curv_init");
  
  // set value
  int icmp = 0;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","x");
  gdcurv->x2d = gdcurv->v3d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","z");
  gdcurv->z2d = gdcurv->v3d + cmp_pos[icmp];
  
  // set pointer
  gdcurv->cmp_pos  = cmp_pos;
  gdcurv->cmp_name = cmp_name;

  return 0;
}

int
grid_init_set(gd_t *gdcurv, char *geometry_file)
{
  FILE *fp = NULL;
  char str[500];
  
  int nx;
  int nz;
  // open geometry file
  if ((fp = fopen(geometry_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open geometry file=%s\n", geometry_file);
     fflush(stdout); exit(1);
  }
  // nx number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nx);
  }
  // nz number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nz);
  }
  
  init_gdcurv(gdcurv,nx,nz);

  // malloc 4 bdry space, read bdry coods
  float *x1;
  float *x2;
  float *z1;
  float *z2;
  x1 = (float *)mem_calloc_1d_float(
            nz*gdcurv->ncmp, 0.0, "bdry_coords");
  x2 = (float *)mem_calloc_1d_float(
            nz*gdcurv->ncmp, 0.0, "bdry_coords");
  z1 = (float *)mem_calloc_1d_float(
            nx*gdcurv->ncmp, 0.0, "bdry_coords");
  z2 = (float *)mem_calloc_1d_float(
            nx*gdcurv->ncmp, 0.0, "bdry_coords");
  // x1 
  for (int k=0; k<nz; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x1+k,(x1+nz)+k);
    }
  }
  // x2 
  for (int k=0; k<nz; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x2+k,(x2+nz)+k);
    }
  }
  // z1 
  for (int i=0; i<nx; i++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",z1+i,(z1+nx)+i);
    }
  }
  // z2 
  for (int i=0; i<nx; i++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",z2+i,(z2+nx)+i);
    }
  }
  // close file and free local pointer
  fclose(fp); 

  check_bdry(x1,x2,z1,z2,nx,nz);

  // read boundry grid coords 
  size_t iptr;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  // x1 i=0
  for (int k=0; k<nz; k++)
  {
    iptr = k*nx;
    x2d[iptr] = x1[k];
    z2d[iptr] = x1[k+nz];
  }
  // x2 i=nx-1
  for (int k=0; k<nz; k++)
  {
    iptr = k*nx + (nx-1);
    x2d[iptr] = x2[k];
    z2d[iptr] = x2[k+nz];
  }
  // z1 k=0
  for (int i=0; i<nx; i++)
  {
    iptr = i;
    x2d[iptr] = z1[i];
    z2d[iptr] = z1[i+nx];
  }
  // z2 k=nz-1
  for (int i=0; i<nx; i++)
  {
    iptr = (nz-1)*nx + i;
    x2d[iptr] = z2[i];
    z2d[iptr] = z2[i+nx];
  }

  free(x1);
  free(x2);
  free(z1);
  free(z2);

  return 0;
}

int
grid_init_set_hyper(gd_t *gdcurv, char *geometry_file, char *step_file)
{
  FILE *fp = NULL;
  char str[500];
  
  int nx;
  int nz;
  int num_step;

  // open step file
  if ((fp = fopen(step_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open step file=%s\n", step_file);
     fflush(stdout); exit(1);
  }
  // number of step
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_step);
  }
  nz = num_step+1; 
  gdcurv->step = (float *)mem_calloc_1d_float(
                          num_step, 0.0, "step length");

  for (int k=0; k<num_step; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f",gdcurv->step + k);
    }
  }
  // close step file and free local pointer
  fclose(fp);

  // open geometry file
  if ((fp = fopen(geometry_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open geometry file=%s\n", geometry_file);
     fflush(stdout); exit(1);
  }
  // nx number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nx);
  }

  init_gdcurv(gdcurv,nx,nz);
 
  size_t iptr;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  for (int i=0; i<nx; i++)
  {
    iptr = i;  // (i,0)
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x2d+iptr,z2d+iptr);
    }
  }
  // close  geometry file and free local pointer
  fclose(fp);

  return 0;
}

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int nx_new = (int) (nx*coef_x);
  int nz_new = (int) (nz*coef_z);
  if(nx_new < nx || nz_new < nz)
  {
    fprintf(stdout,"only support up sample, \
                    nx_new(nz_new) must >= nx(nz)\n");
    exit(1);
  }

  init_gdcurv(gdcurv_new, nx_new, nz_new);
    
  sample_interp(gdcurv_new, gdcurv); 

  return 0;
}

int 
check_bdry(float *x1, float *x2, float *z1, float *z2, int nx, int nz)
{ 
  int ierr = 0;
  float p1_x, p1_z, p2_x, p2_z;
  //  check four corner points
  //  (0,0)
  p1_x = x1[0];  
  p1_z = x1[0+nz];  
  p2_x = z1[0];
  p2_z = z1[0+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (0,0) error, please check x1 and z1 boundary");
    exit(1);
  }

  //  (nx-1,0)
  p1_x = x2[0];  
  p1_z = x2[0+nz];  
  p2_x = z1[nx-1];
  p2_z = z1[nx-1+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (nx-1,0) error, please check x2 and z1 boundary");
    exit(1);
  }

  //  (nx-1,nz-1)
  p1_x = x2[nz-1];  
  p1_z = x2[nz-1+nz];  
  p2_x = z2[nx-1];
  p2_z = z2[nx-1+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (nx-1,nz-1) error, please check x2 and z2 boundary");
    exit(1);
  }

  //  (0,nz-1)
  p1_x = x1[nz-1];  
  p1_z = x1[nz-1+nz];  
  p2_x = z2[0];
  p2_z = z2[0+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (0,nz-1) error, please check  x1 and z2  boundary");
    exit(1);
  }

  return 0;
}

// 2D array flip z direction.  nz-1->0 0->nz-1 i->(nz-1)-i 
int
flip_coord(float *coord, int nx, int nz)
{
  size_t iptr,iptr1;
  float *tmp_coord = NULL;
  tmp_coord = (float *) malloc(nx*nz*sizeof(float));
  // copy data
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++) 
    {
      iptr = k*nx + i;
      tmp_coord[iptr] = coord[iptr];
    }
  }
  // flip coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++) 
    {
      iptr = k*nx + i;
      iptr1 = (nz-1-k)*nx + i;
      coord[iptr] = tmp_coord[iptr1];
    }
  }

  free(tmp_coord);

  return 0;
}

// 2D array permute. transposition (nx,nz) -> (nz,nx)
int
permute_coord(gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  size_t iptr,iptr1;

  float *tmp_coord_x = NULL;
  float *tmp_coord_z = NULL;
  tmp_coord_x = (float *) malloc(nx*nz*sizeof(float));
  tmp_coord_z = (float *) malloc(nx*nz*sizeof(float));
  // copy x, z
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++) 
    {
      iptr = k*nx + i;
      tmp_coord_x[iptr] = x2d[iptr];
      tmp_coord_z[iptr] = z2d[iptr];
    }
  }
  // permute coord, x to z
  // z to x
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++) 
    {
      iptr = i*nz + k;
      iptr1 = k*nx + i;
      z2d[iptr] = tmp_coord_x[iptr1];
      x2d[iptr] = tmp_coord_z[iptr1];
    }
  }

  // NOTE:  
  gdcurv->nx = nz;
  gdcurv->nz = nx;

  free(tmp_coord_x); 
  free(tmp_coord_z); 

  return 0;
}
