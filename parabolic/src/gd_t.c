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
init_gdcurv(gd_t *gdcurv)
{
  //3 dimension, x y and z
  gdcurv->ncmp = 3; 
  // x dimention varies first
  gdcurv->siz_iy = gdcurv->nx; 
  gdcurv->siz_iz = gdcurv->nx*gdcurv->ny; 
  gdcurv->siz_icmp = gdcurv->siz_iz*gdcurv->nz;

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
gd_info_set(gd_t *gdcurv, mympi_t *mympi,
            par_t *par)
{
  int ierr = 0;
  
  if(par->dire_itype == Z_DIRE)
  {
    gdcurv->total_nx = par->number_of_grid_points_x;
    gdcurv->total_ny = par->number_of_grid_points_y;
    gdcurv->total_nz = par->number_of_grid_points_z;
  }
  if(par->dire_itype == Y_DIRE)
  {
    gdcurv->total_nx = par->number_of_grid_points_x;
    gdcurv->total_ny = par->number_of_grid_points_z;
    gdcurv->total_nz = par->number_of_grid_points_y;
  }
  if(par->dire_itype == X_DIRE)
  {
    gdcurv->total_nx = par->number_of_grid_points_z;
    gdcurv->total_ny = par->number_of_grid_points_y;
    gdcurv->total_nz = par->number_of_grid_points_x;
  }

  int number_of_grid_points_x = gdcurv->total_nx;
  int number_of_grid_points_y = gdcurv->total_ny;
  int number_of_grid_points_z = gdcurv->total_nz;

  // determine nj
  int ny_et = number_of_grid_points_y-2;

  int ny_avg  = ny_et / mympi->nproc;
  int ny_left = ny_et % mympi->nproc;

  int nj = ny_avg;

  // not equal divided points given to first ny_left procs
  if (mympi->topoid[0] < ny_left) {
    nj++;
  }
  // global index
  if (mympi->topoid[0]==0) {
    gdcurv->gnj1 = 0;
  } else {
    gdcurv->gnj1 = mympi->topoid[0] * ny_avg;
  }
  if (ny_left != 0) {
    gdcurv->gnj1 += (mympi->topoid[0] < ny_left)? mympi->topoid[0] : ny_left;
  }

  // determine nj
  int ni = number_of_grid_points_x-2;
  int nk = number_of_grid_points_z-2;

  // add ghost point
  int nx = ni + 2;
  int ny = nj + 2;
  int nz = nk + 2;

  gdcurv->ni = ni;
  gdcurv->nj = nj;
  gdcurv->nk = nk;

  gdcurv->nx = nx;
  gdcurv->ny = ny;
  gdcurv->nz = nz;

  gdcurv->ni1 = 1;
  gdcurv->ni2 = gdcurv->ni1 + ni - 1;

  gdcurv->nj1 = 1;
  gdcurv->nj2 = gdcurv->nj1 + nj - 1;

  gdcurv->nk1 = 1;
  gdcurv->nk2 = gdcurv->nk1 + nk - 1;

  // global index end
  gdcurv->gni1 = 0;
  gdcurv->gnk1 = 0;

  gdcurv->gni2 = gdcurv->gni1 + gdcurv->ni - 1;
  gdcurv->gnj2 = gdcurv->gnj1 + gdcurv->nj - 1;
  gdcurv->gnk2 = gdcurv->gnk1 + gdcurv->nk - 1;
  
  return ierr;
}

int
gd_info_print(gd_t *gdcurv, mympi_t *mympi)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "my rank id is %d\n", mympi->myid);
  fprintf(stdout, " nx    = %-10d\n", gdcurv->nx);
  fprintf(stdout, " ny    = %-10d\n", gdcurv->ny);
  fprintf(stdout, " nz    = %-10d\n", gdcurv->nz);
  fprintf(stdout, " ni    = %-10d\n", gdcurv->ni);
  fprintf(stdout, " nj    = %-10d\n", gdcurv->nj);
  fprintf(stdout, " nk    = %-10d\n", gdcurv->nk);

  fprintf(stdout, " ni1   = %-10d\n", gdcurv->ni1);
  fprintf(stdout, " ni2   = %-10d\n", gdcurv->ni2);
  fprintf(stdout, " nj1   = %-10d\n", gdcurv->nj1);
  fprintf(stdout, " nj2   = %-10d\n", gdcurv->nj2);
  fprintf(stdout, " nk1   = %-10d\n", gdcurv->nk1);
  fprintf(stdout, " nk2   = %-10d\n", gdcurv->nk2);

  fprintf(stdout, " ni1_to_glob_phys0   = %-10d\n", gdcurv->gni1);
  fprintf(stdout, " ni2_to_glob_phys0   = %-10d\n", gdcurv->gni2);
  fprintf(stdout, " nj1_to_glob_phys0   = %-10d\n", gdcurv->gnj1);
  fprintf(stdout, " nj2_to_glob_phys0   = %-10d\n", gdcurv->gnj2);
  fprintf(stdout, " nk1_to_glob_phys0   = %-10d\n", gdcurv->gnk1);
  fprintf(stdout, " nk2_to_glob_phys0   = %-10d\n", gdcurv->gnk2);

  fflush(stdout);

  return 0;
}

int
set_output_dir(gd_t *gdcurv, mympi_t *mympi,
               par_t *par)
{
  // output file name
  if(par->dire_itype == Z_DIRE)
  {
    sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", 0, mympi->topoid[0], 0);
  }
  if(par->dire_itype == Y_DIRE)
  {
    sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", 0, 0, mympi->topoid[0]);
  }
  if(par->dire_itype == X_DIRE)
  {
    sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", 0, mympi->topoid[0], 0);
  }

  // output
  sprintf(gdcurv->output_dir, "%s", par->grid_export_dir);

  return 0;
}

int
read_bdry_file(gd_t *gdcurv, par_t *par)
{
  FILE *fp = NULL;
  char str[500];
  char *geometry_file = par->geometry_input_file;
  char *step_file = par->step_input_file;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int gnj1 = gdcurv->gnj1;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  
  int num_step;
  int num_point_x;
  int num_point_y;
  // open step file
  if ((fp = fopen(step_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open step file=%s\n", step_file);
     fflush(stdout); exit(1);
  }
  // number of step
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_step);
  }

  if(nz-num_step != 1) 
  {
    fprintf(stdout,"num_step is error, please check input file");
    fflush(stdout); exit(1);
  }

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
  // nx ny number
  if(par->dire_itype == Z_DIRE || par->dire_itype == Y_DIRE)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d %d",&num_point_x, &num_point_y);
    }
  }
  if(par->dire_itype == X_DIRE)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%d %d",&num_point_y, &num_point_x);
    }
  }

  size_t iptr, iptr1;
  size_t siz_bz = gdcurv->total_nx*gdcurv->total_ny;
  float *bz1 = (float *)mem_calloc_1d_float(siz_bz*3, 0.0, "bottom bdry");
  float *bz2 = (float *)mem_calloc_1d_float(siz_bz*3, 0.0, "top bdry");

  if(par->dire_itype == Z_DIRE)
  {
    // bz1 
    for (int j=0; j<num_point_y; j++) {
      for (int i=0; i<num_point_x; i++) {
        iptr = j*num_point_x+i;  // (i,j)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",bz1+iptr,bz1+siz_bz+iptr,bz1+2*siz_bz+iptr);
        }
      }
    }
    // bz2 
    for (int j=0; j<num_point_y; j++) {
      for (int i=0; i<num_point_x; i++) {
        iptr = j*num_point_x+i;  // (i,j)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",bz2+iptr,bz2+siz_bz+iptr,bz2+2*siz_bz+iptr);
        }
      }
    }
  }

  if(par->dire_itype == Y_DIRE)
  {
    // bz1 
    for (int j=0; j<num_point_y; j++) {
      for (int i=0; i<num_point_x; i++) {
        iptr = j*num_point_x+i;  // (i,j)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",bz1+iptr,bz1+2*siz_bz+iptr,bz1+siz_bz+iptr);
        }
      }
    }
    // bz2 
    for (int j=0; j<num_point_y; j++) {
      for (int i=0; i<num_point_x; i++) {
        iptr = j*num_point_x+i;  // (i,j)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",bz2+iptr,bz2+2*siz_bz+iptr,bz2+siz_bz+iptr);
        }
      }
    }
  }

  if(par->dire_itype == X_DIRE)
  {
    // bz1 
    for (int i=0; i<num_point_x; i++) {
      for (int j=0; j<num_point_y; j++) {
        iptr = i*num_point_y+j;  // (i,j)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",bz1+2*siz_bz+iptr,bz1+siz_bz+iptr,bz1+iptr);
        }
      }
    }
    // bz2 
    for (int i=0; i<num_point_x; i++) {
      for (int j=0; j<num_point_y; j++) {
        iptr = i*num_point_y+j;  // (i,j)
        if (!io_get_nextline(fp,str,500)) {
          sscanf(str,"%f %f %f",bz2+2*siz_bz+iptr,bz2+siz_bz+iptr,bz2+iptr);
        }
      }
    }
  }
  // close geometry file and free local pointer
  fclose(fp);

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  if(par->dire_itype == Z_DIRE || par->dire_itype == Y_DIRE)
  {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        iptr = j*siz_iy+i;  // (i,j,0)
        iptr1 = (gnj1+j)*num_point_x+i;  // (i,j)
        x3d[iptr] = bz1[iptr1];
        y3d[iptr] = bz1[iptr1+siz_bz];
        z3d[iptr] = bz1[iptr1+2*siz_bz];
      }
    }
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        iptr = (nz-1)*siz_iz+j*siz_iy+i;  // (i,j,nz-1)
        iptr1 = (gnj1+j)*num_point_x+i;  // (i,j)
        x3d[iptr] = bz2[iptr1];
        y3d[iptr] = bz2[iptr1+siz_bz];
        z3d[iptr] = bz2[iptr1+2*siz_bz];
      }
    }
  }

  if(par->dire_itype == X_DIRE)
  {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        iptr = j*siz_iy+i;  // (i,j,0)
        iptr1 = i*num_point_y +  (gnj1+j);  // (i,j)
        x3d[iptr] = bz1[iptr1];
        y3d[iptr] = bz1[iptr1+siz_bz];
        z3d[iptr] = bz1[iptr1+2*siz_bz];
      }
    }
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        iptr = (nz-1)*siz_iz+j*siz_iy+i;  // (i,j,nz-1)
        iptr1 = i*num_point_y +  (gnj1+j);  // (i,j)
        x3d[iptr] = bz2[iptr1];
        y3d[iptr] = bz2[iptr1+siz_bz];
        z3d[iptr] = bz2[iptr1+2*siz_bz];
      }
    }
  }

  free(bz1);
  free(bz2);

  return 0;
}


// 3D array permute. transposition (nx,ny,nz) -> (nz,ny,nx)
int
permute_coord_x(gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;

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
        x3d[iptr] = tmp_coord_z[iptr1];
        y3d[iptr] = tmp_coord_y[iptr1];
        z3d[iptr] = tmp_coord_x[iptr1];
      }
    }
  }

  // NOTE: modify nx nz ...  
  gdcurv->nx = nz;
  gdcurv->nz = nx;
  gdcurv->siz_iy = nz;
  gdcurv->siz_iz = nz*ny;
  gdcurv->ni = nk;
  gdcurv->nk = ni;
  gdcurv->ni1 = nk1;
  gdcurv->ni2 = nk2;
  gdcurv->nk1 = ni1;
  gdcurv->nk2 = ni2;
  gdcurv->gni1 = gnk1;
  gdcurv->gnk1 = gni1;

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
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int gnj1 = gdcurv->gnj1;
  int gnk1 = gdcurv->gnk1;

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
  gdcurv->nj = nk;
  gdcurv->nk = nj;
  gdcurv->nj1 = nk1;
  gdcurv->nj2 = nk2;
  gdcurv->nk1 = nj1;
  gdcurv->nk2 = nj2;
  gdcurv->gnj1 = gnk1;
  gdcurv->gnk1 = gnj1;

  free(tmp_coord_x); 
  free(tmp_coord_y); 
  free(tmp_coord_z); 

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
  float *tmp_coord_x = (float *) malloc(siz_icmp*sizeof(float));
  float *tmp_coord_y = (float *) malloc(siz_icmp*sizeof(float));
  float *tmp_coord_z = (float *) malloc(siz_icmp*sizeof(float));
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

int
grid_mesg_init(mympi_t *mympi, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;

  // 3:x,y,z 2: k layer and k+1 layer
  mympi->siz_sbuff_y1 = 3*2*nx;
  mympi->siz_sbuff_y2 = 3*2*nx;

  mympi->siz_rbuff_y1 = 3*2*nx;
  mympi->siz_rbuff_y2 = 3*2*nx;
  

  mympi->siz_sbuff = mympi->siz_sbuff_y1 + mympi->siz_sbuff_y2;

  mympi->siz_rbuff = mympi->siz_rbuff_y1 + mympi->siz_rbuff_y2;

  mympi->sbuff = (float *) malloc(mympi->siz_sbuff*sizeof(MPI_FLOAT));
  mympi->rbuff = (float *) malloc(mympi->siz_rbuff*sizeof(MPI_FLOAT));
  
  int tag[2] = {11, 12};

  mympi->s_reqs = (MPI_Request *) malloc(2*sizeof(MPI_Request));
  mympi->r_reqs = (MPI_Request *) malloc(2*sizeof(MPI_Request));
  // send
  float *sbuff_y1 = mympi->sbuff;
  float *sbuff_y2 = sbuff_y1 + mympi->siz_sbuff_y1;
  MPI_Send_init(sbuff_y1, mympi->siz_sbuff_y1, MPI_FLOAT, mympi->neighid[0], tag[0],
                mympi->topocomm, &(mympi->s_reqs[0]));
  MPI_Send_init(sbuff_y2, mympi->siz_sbuff_y2, MPI_FLOAT, mympi->neighid[1], tag[1],
                mympi->topocomm, &(mympi->s_reqs[1]));

  // recv
  float *rbuff_y1 = mympi->rbuff;
  float *rbuff_y2 = rbuff_y1 + mympi->siz_rbuff_y1;
  MPI_Recv_init(rbuff_y1, mympi->siz_rbuff_y1, MPI_FLOAT, mympi->neighid[0], tag[1],
                mympi->topocomm, &(mympi->r_reqs[0]));
  MPI_Recv_init(rbuff_y2, mympi->siz_rbuff_y2, MPI_FLOAT, mympi->neighid[1], tag[0],
                mympi->topocomm, &(mympi->r_reqs[1]));

  return 0;
}

int
grid_pack_mesg(mympi_t *mympi, gd_t *gdcurv, int k)
{
  int nx = gdcurv->nx;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t iptr, iptr1, iptr2;

  float *sbuff_y1 = mympi->sbuff;
  float *sbuff_y2 = sbuff_y1 + mympi->siz_sbuff_y1;

  // y1 
  for(int i=0; i<nx; i++)
  {
    iptr =  i;
    iptr1 = k*siz_iz + nj1*siz_iy + i;
    iptr2 = (k+1)*siz_iz + nj1*siz_iy + i;
    sbuff_y1[iptr+0*nx] = x3d[iptr1];
    sbuff_y1[iptr+1*nx] = y3d[iptr1];
    sbuff_y1[iptr+2*nx] = z3d[iptr1];
    sbuff_y1[iptr+3*nx] = x3d[iptr2];
    sbuff_y1[iptr+4*nx] = y3d[iptr2];
    sbuff_y1[iptr+5*nx] = z3d[iptr2];
  }

  // y2 
  for(int i=0; i<nx; i++)
  {
    iptr =  i;
    iptr1 = k*siz_iz + nj2*siz_iy + i;
    iptr2 = (k+1)*siz_iz + nj2*siz_iy + i;
    sbuff_y2[iptr+0*nx] = x3d[iptr1];
    sbuff_y2[iptr+1*nx] = y3d[iptr1];
    sbuff_y2[iptr+2*nx] = z3d[iptr1];
    sbuff_y2[iptr+3*nx] = x3d[iptr2];
    sbuff_y2[iptr+4*nx] = y3d[iptr2];
    sbuff_y2[iptr+5*nx] = z3d[iptr2];
  }

  return 0;
}

int
grid_unpack_mesg(mympi_t *mympi, gd_t *gdcurv, int k)
{
  int nx = gdcurv->nx;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int *neighid = mympi->neighid;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t iptr, iptr1, iptr2;

  float *rbuff_y1 = mympi->rbuff;
  float *rbuff_y2 = rbuff_y1 + mympi->siz_rbuff_y1;

  // y1 
  if(neighid[0] != MPI_PROC_NULL)
  {
    for(int i=0; i<nx; i++)
    {
      iptr = i;
      iptr1 = k*siz_iz + (nj1-1)*siz_iy + i;
      iptr2 = (k+1)*siz_iz + (nj1-1)*siz_iy + i;
      x3d[iptr1] = rbuff_y1[iptr+0*nx];
      y3d[iptr1] = rbuff_y1[iptr+1*nx];
      z3d[iptr1] = rbuff_y1[iptr+2*nx];
      x3d[iptr2] = rbuff_y1[iptr+3*nx];
      y3d[iptr2] = rbuff_y1[iptr+4*nx];
      z3d[iptr2] = rbuff_y1[iptr+5*nx];
    }
  }

  // y2
  if(neighid[1] != MPI_PROC_NULL)
  {
    for(int i=0; i<nx; i++)
    {
      iptr = i;
      iptr1 = k*siz_iz + (nj2+1)*siz_iy + i;
      iptr2 = (k+1)*siz_iz + (nj2+1)*siz_iy + i;
      x3d[iptr1] = rbuff_y2[iptr+0*nx];
      y3d[iptr1] = rbuff_y2[iptr+1*nx];
      z3d[iptr1] = rbuff_y2[iptr+2*nx];
      x3d[iptr2] = rbuff_y2[iptr+3*nx];
      y3d[iptr2] = rbuff_y2[iptr+4*nx];
      z3d[iptr2] = rbuff_y2[iptr+5*nx];
    }
  }

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
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;

  for (int k = nk1; k < nk2; k++)
  {
    for (int j = nj1; j < nj2; j++)
    {
      for (int i = ni1; i < ni2; i++)
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
