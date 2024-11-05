#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib_mem.h"
#include "lib_math.h"
#include "gd_t.h"
#include "constants.h"
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
  // NOTE: must use gdcurv->siz_iz
  // because nx*ny*nz is int type
  // when nx,ny,nz is big, nx*ny*nz > int range
  gdcurv->siz_icmp = gdcurv->siz_iz*nz; 
  // malloc grid space 
  gdcurv->v4d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }

  // set value
  int icmp = 0;
  gdcurv->x3d = gdcurv->v4d + icmp * gdcurv->siz_icmp;

  icmp += 1;
  gdcurv->y3d = gdcurv->v4d + icmp * gdcurv->siz_icmp;
  
  icmp += 1;
  gdcurv->z3d = gdcurv->v4d + icmp * gdcurv->siz_icmp;

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

  float *tmp_coord_x = (float *) malloc(siz_icmp*sizeof(float));
  float *tmp_coord_y = (float *) malloc(siz_icmp*sizeof(float));
  float *tmp_coord_z = (float *) malloc(siz_icmp*sizeof(float));
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

  float *tmp_coord_x = (float *) malloc(siz_icmp*sizeof(float));
  float *tmp_coord_y = (float *) malloc(siz_icmp*sizeof(float));
  float *tmp_coord_z = (float *) malloc(siz_icmp*sizeof(float));
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

int
cal_min_dist(gd_t *gdcurv, int *indx_i, int *indx_j, int *indx_k, float *dL_min)
{
  float dL_min_local = 1e10;
  float dL_min_global = 1e10;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;

  for (int k = 1; k < nz-1; k++)
  {
    for (int j = 1; j < ny-1; j++)
    {
      for (int i = 1; i < nx-1; i++)
      {
        size_t iptr = i + j * siz_iy + k * siz_iz;
        float p0[] = { x3d[iptr], y3d[iptr], z3d[iptr] };

        // min L to 8 adjacent planes
        for (int kk = -1; kk <=1; kk = kk+2)
        {
          for (int jj = -1; jj <= 1; jj = jj+2)
          {
            for (int ii = -1; ii <= 1; ii = ii+2) 
            {
              float p1[] = { x3d[iptr-ii], y3d[iptr-ii], z3d[iptr-ii] };
              float p2[] = { x3d[iptr-jj*siz_iy],
                             y3d[iptr-jj*siz_iy],
                             z3d[iptr-jj*siz_iy] };
              float p3[] = { x3d[iptr-kk*siz_iz],
                             y3d[iptr-kk*siz_iz],
                             z3d[iptr-kk*siz_iz] };

              float L = dist_point2plane(p0, p1, p2, p3);

              if (dL_min_local > L) dL_min_local = L;
            }
          }
        }

        if (dL_min_global > dL_min_local) 
        {
          dL_min_global = dL_min_local;
          *dL_min = dL_min_global;
          *indx_i = i;
          *indx_j = j;
          *indx_k = k;
        }
      }
    }
  }

  return 0;
}
