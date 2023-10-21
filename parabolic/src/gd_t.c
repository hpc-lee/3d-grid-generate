#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib_mem.h"
#include "gd_t.h"
#include "constants.h"
#include "io_funcs.h"

int
init_gdcurv(gd_t *gdcurv)
{
  //3 dimension, x y and z
  gdcurv->ncmp = CONST_NDIM; 
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
            par_t *par, int verbose)
{
  int ierr = 0;

  if(par->dire_itype == Z_DIRE)
  {
    gdcurv->total_nx = par->number_of_grid_points_x;
    gdcurv->total_ny = par->number_of_grid_points_y;
    gdcurv->total_nz = par->number_of_grid_points_z;
  }
  // trans y to z, z to y
  if(par->dire_itype == Y_DIRE)
  {
    gdcurv->total_nx = par->number_of_grid_points_x;
    gdcurv->total_ny = par->number_of_grid_points_z;
    gdcurv->total_nz = par->number_of_grid_points_y;
  }
  // trans x to z, z to x
  if(par->dire_itype == X_DIRE)
  {
    gdcurv->total_nx = par->number_of_grid_points_z;
    gdcurv->total_ny = par->number_of_grid_points_y;
    gdcurv->total_nz = par->number_of_grid_points_x;
  }

  int number_of_grid_points_x = gdcurv->total_nx;
  int number_of_grid_points_y = gdcurv->total_ny;
  int number_of_grid_points_z = gdcurv->total_nz;

  // bbry point is ghost point
  // x direction only 1 procs
  // determine ni
  int nx_et = number_of_grid_points_x-2;

  int nx_avg  = nx_et / mympi->nprocx;
  int nx_left = nx_et % mympi->nprocx;

  int ni = nx_avg;

  // not equal divided points given to first ny_left procs
  if (mympi->topoid[0] < nx_left) {
    ni++;
  }
  // global index
  if (mympi->topoid[0]==0) {
    gdcurv->gni1 = 0;
  } else {
    gdcurv->gni1 = mympi->topoid[0] * nx_avg;
  }
  if (nx_left != 0) {
    gdcurv->gni1 += (mympi->topoid[0] < nx_left)? mympi->topoid[0] : nx_left;
  }

  // determine nj
  int ny_et = number_of_grid_points_y-2;

  int ny_avg  = ny_et / mympi->nprocy;
  int ny_left = ny_et % mympi->nprocy;

  int nj = ny_avg;

  // not equal divided points given to first ny_left procs
  if (mympi->topoid[1] < ny_left) {
    nj++;
  }
  // global index
  if (mympi->topoid[1]==0) {
    gdcurv->gnj1 = 0;
  } else {
    gdcurv->gnj1 = mympi->topoid[1] * ny_avg;
  }
  if (ny_left != 0) {
    gdcurv->gnj1 += (mympi->topoid[1] < ny_left)? mympi->topoid[1] : ny_left;
  }

  // determine nk
  int nz_et = number_of_grid_points_z-2;

  int nz_avg  = nz_et / mympi->nprocz;
  int nz_left = nz_et % mympi->nprocz;

  int nk = nz_avg;

  // not equal divided points given to first ny_left procs
  if (mympi->topoid[2] < nz_left) {
    nk++;
  }
  // global index
  if (mympi->topoid[2]==0) {
    gdcurv->gnk1 = 0;
  } else {
    gdcurv->gnk1 = mympi->topoid[2] * nz_avg;
  }
  if (nz_left != 0) {
    gdcurv->gnk1 += (mympi->topoid[2] < nz_left)? mympi->topoid[2] : nz_left;
  }


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
  gdcurv->gni2 = gdcurv->gni1 + gdcurv->ni - 1;
  gdcurv->gnj2 = gdcurv->gnj1 + gdcurv->nj - 1;
  gdcurv->gnk2 = gdcurv->gnk1 + gdcurv->nk - 1;
  
  // x dimention varies first
  gdcurv->siz_iy = nx; 
  gdcurv->siz_iz = nx*ny; 
  gdcurv->siz_icmp = gdcurv->siz_iz*nz;

  return ierr;
}

int
gd_info_print(gd_t *gdcurv, mympi_t *mympi)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "my rank id is %d\n", mympi->myid);
  fprintf(stdout, "topo id is (%d %d %d)\n", 
                   mympi->topoid[0], mympi->topoid[1], mympi->topoid[2]);
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
               par_t *par, int verbose)
{
  // output file name
  if(par->dire_itype == Z_DIRE)
  {
    sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", mympi->topoid[0],mympi->topoid[1],mympi->topoid[2]);
  }
  if(par->dire_itype == Y_DIRE)
  {
    sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", mympi->topoid[0],mympi->topoid[2],mympi->topoid[1]);
  }
  if(par->dire_itype == X_DIRE)
  {
    sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", mympi->topoid[0],mympi->topoid[1],mympi->topoid[2]);
  }

  // output
  sprintf(gdcurv->output_dir, "%s", par->grid_export_dir);

  return 0;
}

int
init_bdry(bdry_t *bdry, par_t *par)
{
  bdry->total_nx = par->number_of_grid_points_x;
  bdry->total_ny = par->number_of_grid_points_y;
  bdry->total_nz = par->number_of_grid_points_z;
  int size_bx = bdry->total_ny*bdry->total_nz;
  int size_by = bdry->total_nx*bdry->total_nz;
  int size_bz = bdry->total_nx*bdry->total_ny;
  int size = 2*(size_bx + size_by + size_bz);
  // malloc bdry space. x y and z coord 
  bdry->var = (float *)mem_calloc_1d_float(size*3, 0.0, "bdry_init");
  if (bdry->var == NULL) {
      fprintf(stderr,"Error: failed to alloc bdry vars\n");
      fflush(stderr);
  }

  bdry->x1 = bdry->var;
  bdry->x2 = bdry->x1+3*size_bx;
  bdry->y1 = bdry->x2+3*size_bx;
  bdry->y2 = bdry->y1+3*size_by;
  bdry->z1 = bdry->y2+3*size_by;
  bdry->z2 = bdry->z1+3*size_bz;

  return 0;
}

int
read_bdry(int myid, bdry_t *bdry, char *geometry_file)
{
  FILE *fp = NULL;
  char str[500];
  
  int total_nx = bdry->total_nx;
  int total_ny = bdry->total_ny;
  int total_nz = bdry->total_nz;

  int size_bx = total_ny*total_nz;
  int size_by = total_nx*total_nz;
  int size_bz = total_nx*total_ny;
  size_t iptr;

  float *x1 = bdry->x1;
  float *x2 = bdry->x2;
  float *y1 = bdry->y1;
  float *y2 = bdry->y2;
  float *z1 = bdry->z1;
  float *z2 = bdry->z2;
  
  int num_point_x;
  int num_point_y;
  int num_point_z;
  // open geometry file
  if ((fp = fopen(geometry_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open geometry file=%s\n", geometry_file);
     fflush(stdout); exit(1);
  }
  // nx number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_point_x);
  }
  // ny number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_point_y);
  }
  // nz number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_point_z);
  }
  // x1 
  for (int k=0; k<total_nz; k++) {
    for (int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",x1+iptr,(x1+size_bx)+iptr,(x1+2*size_bx)+iptr);
      }
    }
  }
  // x2 
  for (int k=0; k<total_nz; k++) {
    for (int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",x2+iptr,(x2+size_bx)+iptr,(x2+2*size_bx)+iptr);
      }
    }
  }
  // y1 
  for (int k=0; k<total_nz; k++) {
    for (int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",y1+iptr,(y1+size_by)+iptr,(y1+2*size_by)+iptr);
      }
    }
  }
  // y2 
  for (int k=0; k<total_nz; k++) {
    for (int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",y2+iptr,(y2+size_by)+iptr,(y2+2*size_by)+iptr);
      }
    }
  }
  // z1 
  for (int j=0; j<total_ny; j++) {
    for (int i=0; i<total_nx; i++)
    {
      iptr = j*total_nx + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",z1+iptr,(z1+size_bz)+iptr,(z1+2*size_bz)+iptr);
      }
    }
  }
  // z2 
  for (int j=0; j<total_ny; j++) {
    for (int i=0; i<total_nx; i++)
    {
      iptr = j*total_nx + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",z2+iptr,(z2+size_bz)+iptr,(z2+2*size_bz)+iptr);
      }
    }
  }
  // close file and free local pointer
  fclose(fp); 
  if(myid == 0)
  {
    check_bdry(x1,x2,y1,y2,z1,z2,total_nx,total_ny,total_nz);
  }

  return 0;
}

// assign 6 bdry point to ghost
int
assign_bdry_coord(gd_t *gdcurv, bdry_t *bdry, mympi_t *mympi)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int gni1 = gdcurv->gni1;
  int gnj1 = gdcurv->gnj1;
  int gnk1 = gdcurv->gnk1;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;
  int gni, gnj, gnk;
  size_t iptr, iptr1;


  float *x1 = bdry->x1;
  float *x2 = bdry->x2;
  float *y1 = bdry->y1;
  float *y2 = bdry->y2;
  float *z1 = bdry->z1;
  float *z2 = bdry->z2;
  size_t size_bx = total_ny*total_nz;
  size_t size_by = total_nx*total_nz;
  size_t size_bz = total_nx*total_ny;

  // bdry x1
  for (int k=0; k<nz; k++)
  {
    for (int j=0; j<ny; j++)
    {
      iptr = k*siz_iz + j*siz_iy;  //(0,j,k) 
      gnk = gnk1 + k;
      gnj = gnj1 + j;
      iptr1 = gnk*total_ny + gnj;

      x3d[iptr] = x1[iptr1+0*size_bx];
      y3d[iptr] = x1[iptr1+1*size_bx];
      z3d[iptr] = x1[iptr1+2*size_bx];
    }
  }
  // bdry x2
  for (int k=0; k<nz; k++)
  {
    for (int j=0; j<ny; j++)
    {
      iptr = k*siz_iz + j*siz_iy + nx-1;  //(nx-1,j,k) 
      gnk = gnk1 + k;
      gnj = gnj1 + j;
      iptr1 = gnk*total_ny + gnj;

      x3d[iptr] = x2[iptr1+0*size_bx];
      y3d[iptr] = x2[iptr1+1*size_bx];
      z3d[iptr] = x2[iptr1+2*size_bx];
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

        x3d[iptr] = y1[iptr1+0*size_by];
        y3d[iptr] = y1[iptr1+1*size_by];
        z3d[iptr] = y1[iptr1+2*size_by];
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

        x3d[iptr] = y2[iptr1+0*size_by];
        y3d[iptr] = y2[iptr1+1*size_by];
        z3d[iptr] = y2[iptr1+2*size_by];
      }
    }
  }
  // bdry z1
  for (int j=0; j<ny; j++)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = j*siz_iy + i;  //(i,j,0) 
      gnj = gnj1 + j; 
      gni = gni1 + i; 
      iptr1 = gnj*total_nx + gni;

      x3d[iptr] = z1[iptr1+0*size_bz];
      y3d[iptr] = z1[iptr1+1*size_bz];
      z3d[iptr] = z1[iptr1+2*size_bz];
    }
  }
  // bdry z2
  for (int j=0; j<ny; j++)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = (nz-1)*siz_iz + j*siz_iy + i;  //(i,j,nz-1) 
      gnj = gnj1 + j; 
      gni = gni1 + i; 
      iptr1 = gnj*total_nx + gni;

      x3d[iptr] = z2[iptr1+0*size_bz];
      y3d[iptr] = z2[iptr1+1*size_bz];
      z3d[iptr] = z2[iptr1+2*size_bz];
    }
  }

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

int
permute_bdry_x(bdry_t *bdry, gd_t *gdcurv)
{
  size_t iptr, iptr1;
  // total_n* has tran in gdcurv
  bdry->total_nx = gdcurv->total_nx;
  bdry->total_ny = gdcurv->total_ny;
  bdry->total_nz = gdcurv->total_nz;

  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;

  float *x1_tmp, *x2_tmp;
  float *z1_tmp, *z2_tmp;
  x1_tmp = bdry->x1;
  x2_tmp = bdry->x2;
  z1_tmp = bdry->z1;
  z2_tmp = bdry->z2;

  bdry->x1 = z1_tmp;
  bdry->x2 = z2_tmp;
  bdry->z1 = x1_tmp;
  bdry->z2 = x2_tmp;

  // copy coord and change sort order
  // old z1(ny*old_nx) -> new x1(new_nz*ny)
  // old_nx = new_nz
  size_t size_bx = total_nz * total_ny;
  float *bx_coord_x = (float *) malloc(sizeof(float)*size_bx);
  float *bx_coord_y = (float *) malloc(sizeof(float)*size_bx);
  float *bx_coord_z = (float *) malloc(sizeof(float)*size_bx);
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr1 = j*total_nz + k;
      iptr  = k*total_ny + j;
      bx_coord_x[iptr] = bdry->x1[iptr1+0*size_bx];
      bx_coord_y[iptr] = bdry->x1[iptr1+1*size_bx];
      bx_coord_z[iptr] = bdry->x1[iptr1+2*size_bx];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      bdry->x1[iptr+0*size_bx] = bx_coord_z[iptr];
      bdry->x1[iptr+1*size_bx] = bx_coord_y[iptr];
      bdry->x1[iptr+2*size_bx] = bx_coord_x[iptr];
    }
  }

  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr1 = j*total_nz + k;
      iptr  = k*total_ny + j;
      bx_coord_x[iptr] = bdry->x2[iptr1+0*size_bx];
      bx_coord_y[iptr] = bdry->x2[iptr1+1*size_bx];
      bx_coord_z[iptr] = bdry->x2[iptr1+2*size_bx];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      bdry->x2[iptr+0*size_bx] = bx_coord_z[iptr];
      bdry->x2[iptr+1*size_bx] = bx_coord_y[iptr];
      bdry->x2[iptr+2*size_bx] = bx_coord_x[iptr];
    }
  }

  // bdry y1 y2
  size_t size_by = total_nz * total_nx;
  float *by_coord_x = (float *) malloc(sizeof(float)*size_by);
  float *by_coord_y = (float *) malloc(sizeof(float)*size_by);
  float *by_coord_z = (float *) malloc(sizeof(float)*size_by);
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = i*total_nz + k;
      iptr  = k*total_nx + i;
      by_coord_x[iptr] = bdry->y1[iptr1+0*size_by];
      by_coord_y[iptr] = bdry->y1[iptr1+1*size_by];
      by_coord_z[iptr] = bdry->y1[iptr1+2*size_by];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      bdry->y1[iptr+0*size_by] = by_coord_z[iptr];
      bdry->y1[iptr+1*size_by] = by_coord_y[iptr];
      bdry->y1[iptr+2*size_by] = by_coord_x[iptr];
    }
  }

  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = i*total_nz + k;
      iptr  = k*total_nx + i;
      by_coord_x[iptr] = bdry->y2[iptr1+0*size_by];
      by_coord_y[iptr] = bdry->y2[iptr1+1*size_by];
      by_coord_z[iptr] = bdry->y2[iptr1+2*size_by];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      bdry->y2[iptr+0*size_by] = by_coord_z[iptr];
      bdry->y2[iptr+1*size_by] = by_coord_y[iptr];
      bdry->y2[iptr+2*size_by] = by_coord_x[iptr];
    }
  }

  // bdry z1 z2
  size_t size_bz = total_ny * total_nx;
  float *bz_coord_x = (float *) malloc(sizeof(float)*size_bz);
  float *bz_coord_y = (float *) malloc(sizeof(float)*size_bz);
  float *bz_coord_z = (float *) malloc(sizeof(float)*size_bz);
  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = i*total_ny + j;
      iptr  = j*total_nx + i;
      bz_coord_x[iptr] = bdry->z1[iptr1+0*size_bz];
      bz_coord_y[iptr] = bdry->z1[iptr1+1*size_bz];
      bz_coord_z[iptr] = bdry->z1[iptr1+2*size_bz];
    }
  }
  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = j*total_nx + i;
      bdry->z1[iptr+0*size_bz] = bz_coord_z[iptr];
      bdry->z1[iptr+1*size_bz] = bz_coord_y[iptr];
      bdry->z1[iptr+2*size_bz] = bz_coord_x[iptr];
    }
  }

  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = i*total_ny + j;
      iptr  = j*total_nx + i;
      bz_coord_x[iptr] = bdry->z2[iptr1+0*size_bz];
      bz_coord_y[iptr] = bdry->z2[iptr1+1*size_bz];
      bz_coord_z[iptr] = bdry->z2[iptr1+2*size_bz];
    }
  }
  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = j*total_nx + i;
      bdry->z2[iptr+0*size_bz] = bz_coord_z[iptr];
      bdry->z2[iptr+1*size_bz] = bz_coord_y[iptr];
      bdry->z2[iptr+2*size_bz] = bz_coord_x[iptr];
    }
  }

  free(bx_coord_x);
  free(bx_coord_y);
  free(bx_coord_z);
  free(by_coord_x);
  free(by_coord_y);
  free(by_coord_z);
  free(bz_coord_x);
  free(bz_coord_y);
  free(bz_coord_z);

  return 0;
}

int
permute_bdry_y(bdry_t *bdry, gd_t *gdcurv)
{
  size_t iptr, iptr1;
  // total_n* has tran in gdcurv
  bdry->total_nx = gdcurv->total_nx;
  bdry->total_ny = gdcurv->total_ny;
  bdry->total_nz = gdcurv->total_nz;

  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;

  float *y1_tmp, *y2_tmp;
  float *z1_tmp, *z2_tmp;
  y1_tmp = bdry->y1;
  y2_tmp = bdry->y2;
  z1_tmp = bdry->z1;
  z2_tmp = bdry->z2;

  bdry->y1 = z1_tmp;
  bdry->y2 = z2_tmp;
  bdry->z1 = y1_tmp;
  bdry->z2 = y2_tmp;

  // copy coord and change sort order
  // old x1(old_nz*old_ny) -> new x1(new_nz*new_ny)
  // old_nx = new_nz old_ny = new_nz
  size_t size_bx = total_nz * total_ny;
  float *bx_coord_x = (float *) malloc(sizeof(float)*size_bx);
  float *bx_coord_y = (float *) malloc(sizeof(float)*size_bx);
  float *bx_coord_z = (float *) malloc(sizeof(float)*size_bx);
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr1 = j*total_nz + k;
      iptr  = k*total_ny + j;
      bx_coord_x[iptr] = bdry->x1[iptr1+0*size_bx];
      bx_coord_y[iptr] = bdry->x1[iptr1+1*size_bx];
      bx_coord_z[iptr] = bdry->x1[iptr1+2*size_bx];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      bdry->x1[iptr+0*size_bx] = bx_coord_x[iptr];
      bdry->x1[iptr+1*size_bx] = bx_coord_z[iptr];
      bdry->x1[iptr+2*size_bx] = bx_coord_y[iptr];
    }
  }

  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr1 = j*total_nz + k;
      iptr  = k*total_ny + j;
      bx_coord_x[iptr] = bdry->x2[iptr1+0*size_bx];
      bx_coord_y[iptr] = bdry->x2[iptr1+1*size_bx];
      bx_coord_z[iptr] = bdry->x2[iptr1+2*size_bx];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int j=0; j<total_ny; j++)
    {
      iptr = k*total_ny + j;
      bdry->x2[iptr+0*size_bx] = bx_coord_x[iptr];
      bdry->x2[iptr+1*size_bx] = bx_coord_z[iptr];
      bdry->x2[iptr+2*size_bx] = bx_coord_y[iptr];
    }
  }

  // bdry y1 y2
  size_t size_by = total_nz * total_nx;
  float *by_coord_x = (float *) malloc(sizeof(float)*size_by);
  float *by_coord_y = (float *) malloc(sizeof(float)*size_by);
  float *by_coord_z = (float *) malloc(sizeof(float)*size_by);
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = k*total_nx + i;
      iptr  = k*total_nx + i;
      by_coord_x[iptr] = bdry->y1[iptr1+0*size_by];
      by_coord_y[iptr] = bdry->y1[iptr1+1*size_by];
      by_coord_z[iptr] = bdry->y1[iptr1+2*size_by];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      bdry->y1[iptr+0*size_by] = by_coord_x[iptr];
      bdry->y1[iptr+1*size_by] = by_coord_z[iptr];
      bdry->y1[iptr+2*size_by] = by_coord_y[iptr];
    }
  }

  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = k*total_nx + i;
      iptr  = k*total_nx + i;
      by_coord_x[iptr] = bdry->y2[iptr1+0*size_by];
      by_coord_y[iptr] = bdry->y2[iptr1+1*size_by];
      by_coord_z[iptr] = bdry->y2[iptr1+2*size_by];
    }
  }
  for(int k=0; k<total_nz; k++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = k*total_nx + i;
      bdry->y2[iptr+0*size_by] = by_coord_x[iptr];
      bdry->y2[iptr+1*size_by] = by_coord_z[iptr];
      bdry->y2[iptr+2*size_by] = by_coord_y[iptr];
    }
  }

  // bdry z1 z2
  size_t size_bz = total_ny * total_nx;
  float *bz_coord_x = (float *) malloc(sizeof(float)*size_bz);
  float *bz_coord_y = (float *) malloc(sizeof(float)*size_bz);
  float *bz_coord_z = (float *) malloc(sizeof(float)*size_bz);
  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = j*total_nx + i;
      iptr  = j*total_nx + i;
      bz_coord_x[iptr] = bdry->z1[iptr1+0*size_bz];
      bz_coord_y[iptr] = bdry->z1[iptr1+1*size_bz];
      bz_coord_z[iptr] = bdry->z1[iptr1+2*size_bz];
    }
  }
  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = j*total_nx + i;
      bdry->z1[iptr+0*size_bz] = bz_coord_x[iptr];
      bdry->z1[iptr+1*size_bz] = bz_coord_z[iptr];
      bdry->z1[iptr+2*size_bz] = bz_coord_y[iptr];
    }
  }

  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr1 = j*total_nx + i;
      iptr  = j*total_nx + i;
      bz_coord_x[iptr] = bdry->z2[iptr1+0*size_bz];
      bz_coord_y[iptr] = bdry->z2[iptr1+1*size_bz];
      bz_coord_z[iptr] = bdry->z2[iptr1+2*size_bz];
    }
  }
  for(int j=0; j<total_ny; j++) {
    for(int i=0; i<total_nx; i++)
    {
      iptr = j*total_nx + i;
      bdry->z2[iptr+0*size_bz] = bz_coord_x[iptr];
      bdry->z2[iptr+1*size_bz] = bz_coord_z[iptr];
      bdry->z2[iptr+2*size_bz] = bz_coord_y[iptr];
    }
  }

  free(bx_coord_x);
  free(bx_coord_y);
  free(bx_coord_z);
  free(by_coord_x);
  free(by_coord_y);
  free(by_coord_z);
  free(bz_coord_x);
  free(bz_coord_y);
  free(bz_coord_z);

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
  int gni1 = gdcurv->gni1;
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
  MPI_Send_init(sbuff_y1, mympi->siz_sbuff_y1, MPI_FLOAT, mympi->neighid[2], tag[0],
                mympi->topocomm, &(mympi->s_reqs[0]));
  MPI_Send_init(sbuff_y2, mympi->siz_sbuff_y2, MPI_FLOAT, mympi->neighid[3], tag[1],
                mympi->topocomm, &(mympi->s_reqs[1]));

  // recv
  float *rbuff_y1 = mympi->rbuff;
  float *rbuff_y2 = rbuff_y1 + mympi->siz_rbuff_y1;
  MPI_Recv_init(rbuff_y1, mympi->siz_rbuff_y1, MPI_FLOAT, mympi->neighid[2], tag[1],
                mympi->topocomm, &(mympi->r_reqs[0]));
  MPI_Recv_init(rbuff_y2, mympi->siz_rbuff_y2, MPI_FLOAT, mympi->neighid[3], tag[0],
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
  if(neighid[2] != MPI_PROC_NULL)
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
  if(neighid[3] != MPI_PROC_NULL)
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
