#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib_mem.h"
#include "gd_t.h"
#include "constants.h"
#include "algebra.h"
#include "io_funcs.h"

int
init_gdcurv(gd_t *gdcurv, int nx, int ny, int nz)
{
  gdcurv->nx = nx;
  gdcurv->ny = ny;
  gdcurv->nz = nz;
  //3 dimension, x y and z
  gdcurv->ncmp = CONST_NDIM; 
  gdcurv->siz_iy   = nx;
  gdcurv->siz_iz   = nx*ny;
  gdcurv->siz_icmp = nx*ny*nz;
  // malloc grid space 
  gdcurv->v4d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v4d == NULL) {
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
  gdcurv->x3d = gdcurv->v4d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","y");
  gdcurv->y3d = gdcurv->v4d + cmp_pos[icmp];
  
  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","z");
  gdcurv->z3d = gdcurv->v4d + cmp_pos[icmp];

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
  int ny;
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
  // ny number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&ny);
  }
  // nz number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nz);
  }
  
  init_gdcurv(gdcurv,nx,ny,nz);

  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;

  // malloc 6 bdry space, read bdry coods
  float *x1;
  float *x2;
  float *y1;
  float *y2;
  float *z1;
  float *z2;
  x1 = (float *)mem_calloc_1d_float(
            nz*ny*gdcurv->ncmp, 0.0, "bdry_coords");
  x2 = (float *)mem_calloc_1d_float(
            nz*ny*gdcurv->ncmp, 0.0, "bdry_coords");
  y1 = (float *)mem_calloc_1d_float(
            nz*nx*gdcurv->ncmp, 0.0, "bdry_coords");
  y2 = (float *)mem_calloc_1d_float(
            nz*nx*gdcurv->ncmp, 0.0, "bdry_coords");
  z1 = (float *)mem_calloc_1d_float(
            ny*nx*gdcurv->ncmp, 0.0, "bdry_coords");
  z2 = (float *)mem_calloc_1d_float(
            ny*nx*gdcurv->ncmp, 0.0, "bdry_coords");

  size_t iptr,size;
  // x1 
  for (int k=0; k<nz; k++) 
  {
    for (int j=0; j<ny; j++)
    {
      iptr = k*ny+j;
      size = ny*nz;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",x1+iptr,(x1+1*size)+iptr,(x1+2*size)+iptr);
      }
    }
  }
  // x2 
  for (int k=0; k<nz; k++) 
  {
    for (int j=0; j<ny; j++)
    {
      iptr = k*ny+j;
      size = ny*nz;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",x2+iptr,(x2+1*size)+iptr,(x2+2*size)+iptr);
      }
    }
  }
  // y1 
  for (int k=0; k<nz; k++) 
  {
    for (int i=0; i<nx; i++)
    {
      iptr = k*nx+i;
      size = nx*nz;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",y1+iptr,(y1+1*size)+iptr,(y1+2*size)+iptr);
      }
    }
  }
  // y2 
  for (int k=0; k<nz; k++) 
  {
    for (int i=0; i<nx; i++)
    {
      iptr = k*nx+i;
      size = nx*nz;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",y2+iptr,(y2+1*size)+iptr,(y2+2*size)+iptr);
      }
    }
  }
  // z1 
  for (int j=0; j<ny; j++) 
  {
    for (int i=0; i<nx; i++)
    {
      iptr = j*nx+i;
      size = nx*ny;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",z1+iptr,(z1+1*size)+iptr,(z1+2*size)+iptr);
      }
    }
  }
  // z2 
  for (int j=0; j<ny; j++) 
  {
    for (int i=0; i<nx; i++)
    {
      iptr = j*nx+i;
      size = nx*ny;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",z2+iptr,(z2+1*size)+iptr,(z2+2*size)+iptr);
      }
    }
  }
  // close file and free local pointer
  fclose(fp); 

  check_bdry(x1,x2,y1,y2,z1,z2,nx,ny,nz);

  // read boundry grid coords 
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t iptr1;
  // x1 i=0
  for (int k=0; k<nz; k++)
  {
    for (int j=0; j<ny; j++)
    {
      iptr = k*siz_iz + j*siz_iy + 0;
      iptr1 = k*ny + j;
      size = ny*nz;
      x3d[iptr] = x1[iptr1];
      y3d[iptr] = x1[iptr1+1*size];
      z3d[iptr] = x1[iptr1+2*size];
    }
  }
  // x2 i=nx-1
  for (int k=0; k<nz; k++)
  {
    for (int j=0; j<ny; j++)
    {
      iptr = k*siz_iz + j*siz_iy + (nx-1);
      iptr1 = k*ny + j;
      size = ny*nz;
      x3d[iptr] = x2[iptr1];
      y3d[iptr] = x2[iptr1+1*size];
      z3d[iptr] = x2[iptr1+2*size];
    }
  }
  // y1 j=0
  for (int k=0; k<nz; k++)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = k*siz_iz + 0*siz_iy + i;
      iptr1 = k*nx + i;
      size = nx*nz;
      x3d[iptr] = y1[iptr1];
      y3d[iptr] = y1[iptr1+1*size];
      z3d[iptr] = y1[iptr1+2*size];
    }
  }
  // y2 j=ny-1
  for (int k=0; k<nz; k++)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = k*siz_iz + (ny-1)*siz_iy + i;
      iptr1 = k*nx + i;
      size = nx*nz;
      x3d[iptr] = y2[iptr1];
      y3d[iptr] = y2[iptr1+1*size];
      z3d[iptr] = y2[iptr1+2*size];
    }
  }
  // z1 k=0
  for (int j=0; j<ny; j++)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = 0*siz_iz + j*siz_iy + i;
      iptr1 = j*nx + i;
      size = nx*ny;
      x3d[iptr] = z1[iptr1];
      y3d[iptr] = z1[iptr1+1*size];
      z3d[iptr] = z1[iptr1+2*size];
    }
  }
  // z2 k=nz-1
  for (int j=0; j<ny; j++)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = (nz-1)*siz_iz + j*siz_iy + i;
      iptr1 = j*nx + i;
      size = nx*ny;
      x3d[iptr] = z2[iptr1];
      y3d[iptr] = z2[iptr1+1*size];
      z3d[iptr] = z2[iptr1+2*size];
    }
  }

  free(x1);
  free(x2);
  free(y1);
  free(y2);
  free(z1);
  free(z2);

  return 0;
}

int
grid_init_set_hyper(gd_t *gdcurv, par_t *par)
{
  FILE *fp = NULL;
  char str[500];
  char *geometry_file = par->geometry_input_file;
  char *step_file = par->step_input_file;
  int dire_itype = par->dire_itype;
  int nx;
  int ny;
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

  if(dire_itype == Z_DIRE) nz = num_step+1;
  if(dire_itype == Y_DIRE) ny = num_step+1;
  if(dire_itype == X_DIRE) nx = num_step+1;
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
  if(dire_itype == Z_DIRE)
  {
    // nx number
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d",&nx);
    }
    // ny number
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d",&ny);
    }

    init_gdcurv(gdcurv,nx,ny,nz);
 
    size_t iptr;
    float *x3d = gdcurv->x3d;
    float *y3d = gdcurv->y3d;
    float *z3d = gdcurv->z3d;
    size_t siz_iy = gdcurv->siz_iy;
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++)
      {
        iptr = j*siz_iy + i;  // (i,j,0)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",x3d+iptr,y3d+iptr,z3d+iptr);
        }
      }
    }
  }
  if(dire_itype == Y_DIRE)
  {
    // nx number
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d",&nx);
    }
    // ny number
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d",&nz);
    }

    init_gdcurv(gdcurv,nx,ny,nz);
 
    size_t iptr;
    float *x3d = gdcurv->x3d;
    float *y3d = gdcurv->y3d;
    float *z3d = gdcurv->z3d;
    size_t siz_iz = gdcurv->siz_iz;
    for (int k=0; k<nz; k++) {
      for (int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + i;  // (i,0,k)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",x3d+iptr,y3d+iptr,z3d+iptr);
        }
      }
    }
    permute_coord_y(gdcurv);
  }
  if(dire_itype == X_DIRE)
  {
    // nx number
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d",&ny);
    }
    // ny number
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d",&nz);
    }

    init_gdcurv(gdcurv,nx,ny,nz);
 
    size_t iptr;
    float *x3d = gdcurv->x3d;
    float *y3d = gdcurv->y3d;
    float *z3d = gdcurv->z3d;
    size_t siz_iy = gdcurv->siz_iy;
    size_t siz_iz = gdcurv->siz_iz;
    for (int k=0; k<nz; k++) {
      for (int j=0; j<ny; j++)
      {
        iptr = k*siz_iz + j*siz_iy;  // (0,j,k)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",x3d+iptr,y3d+iptr,z3d+iptr);
        }
      }
    }
    permute_coord_x(gdcurv);
  }
  // close  geometry file and free local pointer
  fclose(fp);

  return 0;
}

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_y, float coef_z)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int nx_new = (int) (nx*coef_x);
  int ny_new = (int) (ny*coef_y);
  int nz_new = (int) (nz*coef_z);
  if(nx_new < nx || ny_new < ny || nz_new < nz)
  {
    fprintf(stdout,"only support up sample, \
                    nx_new(ny_new,nz_new) must >= nx(ny,nz)\n");
    exit(1);
  }

  init_gdcurv(gdcurv_new, nx_new, ny_new, nz_new);
    
  sample_interp(gdcurv_new, gdcurv); 

  return 0;
}

int 
check_bdry(float *x1, float *x2, float *y1, float *y2, float *z1, float *z2,
           int nx, int ny, int nz)
{ 
  int ierr = 0;
  size_t iptr1, iptr2, size1, size2;
  float dif_x,dif_y,dif_z,dif;
  // check 12 edges line
  // 1 check bdry y1 z1
  for(int i=0; i<nx; i++)
  {
    iptr1 = 0*nx+i; // y1
    iptr2 = 0*nx+i; // z1
    size1 = nz*nx;  // y1
    size2 = ny*nx;  // z1
    dif_x = y1[iptr1] - z1[iptr2];
    dif_y = y1[iptr1+1*size1] - z1[iptr2+1*size2];
    dif_z = y1[iptr1+2*size1] - z1[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
  
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge y1 z1, please check y1 and z1 boundary\n");
      fprintf(stdout, "point %d is error\n",i+1);
      exit(1);
    }
  }
  // 2 check bdry y1 z2
  for(int i=0; i<nx; i++)
  {
    iptr1 = (nz-1)*nx+i; // y1
    iptr2 = 0*nx+i; // z2
    size1 = nz*nx;  // y1
    size2 = ny*nx;  // z2
    dif_x = y1[iptr1] - z2[iptr2];
    dif_y = y1[iptr1+1*size1] - z2[iptr2+1*size2];
    dif_z = y1[iptr1+2*size1] - z2[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge y1 z2, please check y1 and z2 boundary\n");
      fprintf(stdout, "point %d is error\n",i+1);
      exit(1);
    }
  }
  // 3 check bdry y2 z1
  for(int i=0; i<nx; i++)
  {
    iptr1 = 0*nx+i; // y2
    iptr2 = (ny-1)*nx+i; // z1
    size1 = nz*nx;  // y2
    size2 = ny*nx;  // z1
    dif_x = y2[iptr1] - z1[iptr2];
    dif_y = y2[iptr1+1*size1] - z1[iptr2+1*size2];
    dif_z = y2[iptr1+2*size1] - z1[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge y2 z1, please check y2 and z1 boundary\n");
      fprintf(stdout, "point %d is error\n",i+1);
      exit(1);
    }
  }
  // 4 check bdry y2 z2
  for(int i=0; i<nx; i++)
  {
    iptr1 = (nz-1)*nx+i; // y2
    iptr2 = (ny-1)*nx+i; // z2
    size1 = nz*nx;  // y2
    size2 = ny*nx;  // z2
    dif_x = y2[iptr1] - z2[iptr2];
    dif_y = y2[iptr1+1*size1] - z2[iptr2+1*size2];
    dif_z = y2[iptr1+2*size1] - z2[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge y2 z2, please check y2 and z2 boundary\n");
      fprintf(stdout, "point %d is error\n",i+1);
      exit(1);
    }
  }

  // 5 check bdry x1 z1
  for(int j=0; j<ny; j++)
  {
    iptr1 = 0*ny+j; // x1
    iptr2 = j*nx+0; // z1
    size1 = nz*ny;  // x1
    size2 = ny*nx;  // z1
    dif_x = x1[iptr1] - z1[iptr2];
    dif_y = x1[iptr1+1*size1] - z1[iptr2+1*size2];
    dif_z = x1[iptr1+2*size1] - z1[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x1 z1, please check x1 and z1 boundary\n");
      fprintf(stdout, "point %d is error\n",j+1);
      exit(1);
    }
  }
  // 6 check bdry x1 z2
  for(int j=0; j<ny; j++)
  {
    iptr1 = (nz-1)*ny+j; // x1
    iptr2 = j*nx+0; // z2
    size1 = nz*ny;  // x1
    size2 = ny*nx;  // z2
    dif_x = x1[iptr1] - z2[iptr2];
    dif_y = x1[iptr1+1*size1] - z2[iptr2+1*size2];
    dif_z = x1[iptr1+2*size1] - z2[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x1 z2, please check x1 and z2 boundary\n");
      fprintf(stdout, "point %d is error\n",j+1);
      exit(1);
    }
  }
  // 7 check bdry x2 z1
  for(int j=0; j<ny; j++)
  {
    iptr1 = 0*ny+j; // x2
    iptr2 = j*nx+(nx-1); // z1
    size1 = nz*ny;  // x2
    size2 = ny*nx;  // z1
    dif_x = x2[iptr1] - z1[iptr2];
    dif_y = x2[iptr1+1*size1] - z1[iptr2+1*size2];
    dif_z = x2[iptr1+2*size1] - z1[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x2 z1, please check x2 and z1 boundary\n");
      fprintf(stdout, "point %d is error\n",j+1);
      exit(1);
    }
  }
  // 8 check bdry x2 z2
  for(int j=0; j<ny; j++)
  {
    iptr1 = (nz-1)*ny+j; // x2
    iptr2 = j*nx+(nx-1); // z2
    size1 = nz*ny;  // x2
    size2 = ny*nx;  // z2
    dif_x = x2[iptr1] - z2[iptr2];
    dif_y = x2[iptr1+1*size1] - z2[iptr2+1*size2];
    dif_z = x2[iptr1+2*size1] - z2[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x2 z2, please check x2 and z2 boundary\n");
      fprintf(stdout, "point %d is error\n",j+1);
      exit(1);
    }
  }

  // 9 check bdry x1 y1
  for(int k=0; k<nz; k++)
  {
    iptr1 = k*ny+0; // x1
    iptr2 = k*nx+0; // y1
    size1 = nz*ny;  // x1
    size2 = nz*nx;  // y1
    dif_x = x1[iptr1] - y1[iptr2];
    dif_y = x1[iptr1+1*size1] - y1[iptr2+1*size2];
    dif_z = x1[iptr1+2*size1] - y1[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x1 y1, please check x1 and y1 boundary\n");
      fprintf(stdout, "point %d is error\n",k+1);
      exit(1);
    }
  }
  // 10 check bdry x1 y2
  for(int k=0; k<nz; k++)
  {
    iptr1 = k*ny+(ny-1); // x1
    iptr2 = k*nx+0; // y2
    size1 = nz*ny;  // x1
    size2 = nz*nx;  // y2
    dif_x = x1[iptr1] - y2[iptr2];
    dif_y = x1[iptr1+1*size1] - y2[iptr2+1*size2];
    dif_z = x1[iptr1+2*size1] - y2[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x1 y2, please check x1 and y2 boundary\n");
      fprintf(stdout, "point %d is error\n",k+1);
      exit(1);
    }
  }
  // 11 check bdry x2 y1
  for(int k=0; k<nz; k++)
  {
    iptr1 = k*ny+0; // x2
    iptr2 = k*nx+(nx-1); // y1
    size1 = nz*ny;  // x2
    size2 = nz*nx;  // y1
    dif_x = x2[iptr1] - y1[iptr2];
    dif_y = x2[iptr1+1*size1] - y1[iptr2+1*size2];
    dif_z = x2[iptr1+2*size1] - y1[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x2 y1, please check x2 and y1 boundary\n");
      fprintf(stdout, "point %d is error\n",k+1);
      exit(1);
    }
  }
  // 12 check bdry x2 y2
  for(int k=0; k<nz; k++)
  {
    iptr1 = k*ny+(ny-1); // x2
    iptr2 = k*nx+(nx-1); // y2
    size1 = nz*ny;  // x2
    size2 = nz*nx;  // y2
    dif_x = x2[iptr1] - y2[iptr2];
    dif_y = x2[iptr1+1*size1] - y2[iptr2+1*size2];
    dif_z = x2[iptr1+2*size1] - y2[iptr2+2*size2];
    dif = fabs(dif_x) + fabs(dif_y) + fabs(dif_z);
    if(dif > 0.00001)
    {
      ierr = 1;
      fprintf(stdout, "edge x2 y2, please check x2 and y2 boundary\n");
      fprintf(stdout, "point %d is error\n",k+1);
      exit(1);
    }
  }
  return 0;
}

// 3D array flip z direction.  nz-1->0 0->nz-1 i->(nz-1)-i 
int
flip_coord_z(gd_t *gdcurv)
{
  size_t iptr,iptr1;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t siz_icmp = gdcurv->siz_icmp;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *tmp_coord_x = NULL;
  float *tmp_coord_y = NULL;
  float *tmp_coord_z = NULL;
  tmp_coord_x = (float *) malloc(siz_icmp*sizeof(float));
  tmp_coord_y = (float *) malloc(siz_icmp*sizeof(float));
  tmp_coord_z = (float *) malloc(siz_icmp*sizeof(float));
  // copy data
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) 
      {
        iptr = k*siz_iz + j*siz_iy + i;
        tmp_coord_x[iptr] = x3d[iptr];
        tmp_coord_y[iptr] = y3d[iptr];
        tmp_coord_z[iptr] = z3d[iptr];
      }
    }
  }
  // flip coord
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) 
      {
        iptr = k*siz_iz + j*siz_iy + i;
        iptr1 = (nz-1-k)*siz_iz + j*siz_iy  + i;
        x3d[iptr] = tmp_coord_x[iptr1];
        y3d[iptr] = tmp_coord_y[iptr1];
        z3d[iptr] = tmp_coord_z[iptr1];
      }
    }
  }

  free(tmp_coord_x);
  free(tmp_coord_y);
  free(tmp_coord_z);

  return 0;
}

// 3D array permute. transposition (nx,ny,nz) -> (nz,ny,nx)
int
permute_coord_x(gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t siz_icmp = gdcurv->siz_icmp;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  size_t iptr,iptr1;

  float *tmp_coord_x = NULL;
  float *tmp_coord_y = NULL;
  float *tmp_coord_z = NULL;
  tmp_coord_x = (float *) malloc(siz_icmp*sizeof(float));
  tmp_coord_y = (float *) malloc(siz_icmp*sizeof(float));
  tmp_coord_z = (float *) malloc(siz_icmp*sizeof(float));
  // copy x, y, z
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) 
      {
        iptr = k*siz_iz + j*siz_iy + i;
        tmp_coord_x[iptr] = x3d[iptr];
        tmp_coord_y[iptr] = y3d[iptr];
        tmp_coord_z[iptr] = z3d[iptr];
      }
    }
  }
  // permute coord, x to z
  // z to x
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) 
      {
        // NOTE: after trans, size_iz = (ny*nz)
        // size_iy = nz
        iptr  = i*(ny*nz) + j*nz + k;  
        iptr1 = k*siz_iz + j*siz_iy + i;
        z3d[iptr] = tmp_coord_x[iptr1];
        x3d[iptr] = tmp_coord_z[iptr1];
        y3d[iptr] = tmp_coord_y[iptr1];
      }
    }
  }

  // NOTE: modify nx nz ...  
  gdcurv->nx = nz;
  gdcurv->nz = nx;
  gdcurv->siz_iy = nz;
  gdcurv->siz_iz = nz*ny;

  free(tmp_coord_x); 
  free(tmp_coord_y); 
  free(tmp_coord_z); 

  return 0;
}
// 3D array permute. transposition (nx,ny,nz) -> (nx,nz,ny)
int
permute_coord_y(gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t siz_icmp = gdcurv->siz_icmp;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  size_t iptr,iptr1;

  float *tmp_coord_x = NULL;
  float *tmp_coord_y = NULL;
  float *tmp_coord_z = NULL;
  tmp_coord_x = (float *) malloc(siz_icmp*sizeof(float));
  tmp_coord_y = (float *) malloc(siz_icmp*sizeof(float));
  tmp_coord_z = (float *) malloc(siz_icmp*sizeof(float));
  // copy x, y, z
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) 
      {
        iptr = k*siz_iz + j*siz_iy + i;
        tmp_coord_x[iptr] = x3d[iptr];
        tmp_coord_y[iptr] = y3d[iptr];
        tmp_coord_z[iptr] = z3d[iptr];
      }
    }
  }
  // permute coord, x to z
  // z to x
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) 
      {
        // NOTE: after trans, size_iz = (nx*nz)
        // size_iy = nx
        iptr  = j*(nx*nz) + k*nx + i;
        iptr1 = k*siz_iz + j*siz_iy + i;
        z3d[iptr] = tmp_coord_y[iptr1];
        y3d[iptr] = tmp_coord_z[iptr1];
        x3d[iptr] = tmp_coord_x[iptr1];
      }
    }
  }

  // NOTE:  
  gdcurv->ny = nz;
  gdcurv->nz = ny;
  gdcurv->siz_iy = nx;
  gdcurv->siz_iz = nx*nz;

  free(tmp_coord_x); 
  free(tmp_coord_y); 
  free(tmp_coord_z); 

  return 0;
}