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
init_gdcurv(gd_t *gdcurv)
{
  // 3 dimension, x y and z
  gdcurv->ncmp = CONST_NDIM; 
  // malloc grid space
  gdcurv->v4d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // set value
  int icmp;
  icmp = 0;
  gdcurv->x3d = gdcurv->v4d + icmp * gdcurv->siz_icmp;

  icmp = 1;
  gdcurv->y3d = gdcurv->v4d + icmp * gdcurv->siz_icmp;

  icmp = 2;
  gdcurv->z3d = gdcurv->v4d + icmp * gdcurv->siz_icmp;

  // malloc temp grid space
  gdcurv->v4d_tmp = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v4d_tmp == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // set value
  icmp = 0;
  gdcurv->x3d_tmp = gdcurv->v4d_tmp + icmp * gdcurv->siz_icmp;

  icmp = 1;
  gdcurv->y3d_tmp = gdcurv->v4d_tmp + icmp * gdcurv->siz_icmp;

  icmp = 2;
  gdcurv->z3d_tmp = gdcurv->v4d_tmp + icmp * gdcurv->siz_icmp;

  return 0;
}

int
init_bdry(bdry_t *bdry, par_t *par)
{
  bdry->total_nx = par->number_of_grid_points_x;
  bdry->total_ny = par->number_of_grid_points_y;
  bdry->total_nz = par->number_of_grid_points_z;
  int siz_bx = bdry->total_ny*bdry->total_nz;
  int siz_by = bdry->total_nx*bdry->total_nz;
  int siz_bz = bdry->total_nx*bdry->total_ny;
  int size = 2*(siz_bx + siz_by + siz_bz);
  // malloc grid space. x y and z coord 
  bdry->var = (float *)mem_calloc_1d_float(size*3, 0.0, "bdry_init");
  if (bdry->var == NULL) {
      fprintf(stderr,"Error: failed to alloc bdry vars\n");
      fflush(stderr);
  }

  bdry->x1 = bdry->var;
  bdry->x2 = bdry->x1+3*siz_bx;
  bdry->y1 = bdry->x2+3*siz_bx;
  bdry->y2 = bdry->y1+3*siz_by;
  bdry->z1 = bdry->y2+3*siz_by;
  bdry->z2 = bdry->z1+3*siz_bz;

  return 0;
}

int
read_bdry(int myid, bdry_t *bdry, char *geometry_file)
{
  FILE *fp = NULL;
  char str[500];

  int siz_bx = bdry->total_ny*bdry->total_nz;
  int siz_by = bdry->total_nx*bdry->total_nz;
  int siz_bz = bdry->total_nx*bdry->total_ny;
  int iptr;

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
  for (int k=0; k<num_point_z; k++) {
    for (int j=0; j<num_point_y; j++)
    {
      iptr = k*num_point_y + j;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",x1+iptr,(x1+siz_bx)+iptr,(x1+2*siz_bx)+iptr);
      }
    }
  }
  // x2 
  for (int k=0; k<num_point_z; k++) {
    for (int j=0; j<num_point_y; j++)
    {
      iptr = k*num_point_y + j;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",x2+iptr,(x2+siz_bx)+iptr,(x2+2*siz_bx)+iptr);
      }
    }
  }
  // y1 
  for (int k=0; k<num_point_z; k++) {
    for (int i=0; i<num_point_x; i++)
    {
      iptr = k*num_point_x + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",y1+iptr,(y1+siz_by)+iptr,(y1+2*siz_by)+iptr);
      }
    }
  }
  // y2 
  for (int k=0; k<num_point_z; k++) {
    for (int i=0; i<num_point_x; i++)
    {
      iptr = k*num_point_x + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",y2+iptr,(y2+siz_by)+iptr,(y2+2*siz_by)+iptr);
      }
    }
  }
  // z1 
  for (int j=0; j<num_point_y; j++) {
    for (int i=0; i<num_point_x; i++)
    {
      iptr = j*num_point_x + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",z1+iptr,(z1+siz_bz)+iptr,(z1+2*siz_bz)+iptr);
      }
    }
  }
  // z2 
  for (int j=0; j<num_point_y; j++) {
    for (int i=0; i<num_point_x; i++)
    {
      iptr = j*num_point_x + i;
      if (!io_get_nextline(fp,str,500)) {
        sscanf(str,"%f %f %f",z2+iptr,(z2+siz_bz)+iptr,(z2+2*siz_bz)+iptr);
      }
    }
  }
  // close file and free local pointer
  fclose(fp); 
  if(myid == 0)
  {
    check_bdry(x1,x2,y1,y2,z1,z2,num_point_x,num_point_y,num_point_z);
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
gd_info_set(gd_t *gdcurv,
            mympi_t *mympi,
            par_t *par)
           
{
  int ierr = 0;

  gdcurv->total_nx = par->number_of_grid_points_x;
  gdcurv->total_ny = par->number_of_grid_points_y;
  gdcurv->total_nz = par->number_of_grid_points_z;

  int number_of_grid_points_x = par->number_of_grid_points_x;
  int number_of_grid_points_y = par->number_of_grid_points_y;
  int number_of_grid_points_z = par->number_of_grid_points_z;

  // 3 point center difference
  // ghost number is 1

  // determine ni
  // bbry point to set ghost
  int nx_et = number_of_grid_points_x-2;

  // partition into average plus left at last
  int nx_avg  = nx_et / mympi->nprocx;
  int nx_left = nx_et % mympi->nprocx;


  // default set to average value
  int ni = nx_avg;

  // first nx_left node add one more point
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

  // not equal divided points given to first nz_left procs
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
  gdcurv->siz_icmp = nx*ny*nz;

  return ierr;
}

int
set_output_dir(gd_t *gdcurv,
               mympi_t *mympi,
               char *output_dir)
{
  // output file name
  sprintf(gdcurv->fname_part,"px%d_py%d_pz%d", mympi->topoid[0],mympi->topoid[1],mympi->topoid[2]);

  // output
  sprintf(gdcurv->output_dir, "%s", output_dir);

  return 0;
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
gd_curv_coord_exchange(gd_t *gdcurv, int *neighid, MPI_Comm topocomm)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  size_t s_iptr;
  size_t r_iptr;

  MPI_Status status;
  MPI_Datatype DTypeXL, DTypeYL, DTypeZL;

  MPI_Type_vector(ny*nz,
                  1,
                  nx,
                  MPI_FLOAT,
                  &DTypeXL);
  MPI_Type_vector(nz,
                  nx,
                  nx*ny,
                  MPI_FLOAT,
                  &DTypeYL);
  MPI_Type_vector(1,
                  nx*ny,
                  1,
                  MPI_FLOAT,
                  &DTypeZL);
  MPI_Type_commit(&DTypeXL);
  MPI_Type_commit(&DTypeYL);
  MPI_Type_commit(&DTypeZL);

  // x1 to x2
  s_iptr = ni1;
  r_iptr = ni2+1;
  MPI_Sendrecv(&x3d[s_iptr],1,DTypeXL,neighid[0],110,
               &x3d[r_iptr],1,DTypeXL,neighid[1],110,
               topocomm,&status);
  MPI_Sendrecv(&y3d[s_iptr],1,DTypeXL,neighid[0],110,
               &y3d[r_iptr],1,DTypeXL,neighid[1],110,
               topocomm,&status);
  MPI_Sendrecv(&z3d[s_iptr],1,DTypeXL,neighid[0],110,
               &z3d[r_iptr],1,DTypeXL,neighid[1],110,
               topocomm,&status);
  // x2 to x1
  s_iptr = ni2;
  r_iptr = ni1-1;
  MPI_Sendrecv(&x3d[s_iptr],1,DTypeXL,neighid[1],120,
               &x3d[r_iptr],1,DTypeXL,neighid[0],120,
               topocomm,&status);
  MPI_Sendrecv(&y3d[s_iptr],1,DTypeXL,neighid[1],120,
               &y3d[r_iptr],1,DTypeXL,neighid[0],120,
               topocomm,&status);
  MPI_Sendrecv(&z3d[s_iptr],1,DTypeXL,neighid[1],120,
               &z3d[r_iptr],1,DTypeXL,neighid[0],120,
               topocomm,&status);
  // y1 to y2
  s_iptr = nj1 * siz_iy;
  r_iptr = (nj2+1) * siz_iy;
  MPI_Sendrecv(&x3d[s_iptr],1,DTypeYL,neighid[2],210,
               &x3d[r_iptr],1,DTypeYL,neighid[3],210,
               topocomm,&status);
  MPI_Sendrecv(&y3d[s_iptr],1,DTypeYL,neighid[2],210,
               &y3d[r_iptr],1,DTypeYL,neighid[3],210,
               topocomm,&status);
  MPI_Sendrecv(&z3d[s_iptr],1,DTypeYL,neighid[2],210,
               &z3d[r_iptr],1,DTypeYL,neighid[3],210,
               topocomm,&status);
  // y2 to y1
  s_iptr = nj2 * siz_iy;
  r_iptr = (nj1-1) * siz_iy;
  MPI_Sendrecv(&x3d[s_iptr],1,DTypeYL,neighid[3],220,
               &x3d[r_iptr],1,DTypeYL,neighid[2],220,
               topocomm,&status);
  MPI_Sendrecv(&y3d[s_iptr],1,DTypeYL,neighid[3],220,
               &y3d[r_iptr],1,DTypeYL,neighid[2],220,
               topocomm,&status);
  MPI_Sendrecv(&z3d[s_iptr],1,DTypeYL,neighid[3],220,
               &z3d[r_iptr],1,DTypeYL,neighid[2],220,
               topocomm,&status);
  // z1 to z2
  s_iptr = nk1 * siz_iz;        
  r_iptr = (nk2+1) * siz_iz;    
  MPI_Sendrecv(&x3d[s_iptr],1,DTypeZL,neighid[4],310,
               &x3d[r_iptr],1,DTypeZL,neighid[5],310,
               topocomm,&status);
  MPI_Sendrecv(&y3d[s_iptr],1,DTypeZL,neighid[4],310,
               &y3d[r_iptr],1,DTypeZL,neighid[5],310,
               topocomm,&status);
  MPI_Sendrecv(&z3d[s_iptr],1,DTypeZL,neighid[4],310,
               &z3d[r_iptr],1,DTypeZL,neighid[5],310,
               topocomm,&status);
  // z2 to z1
  s_iptr = nk2 * siz_iz;
  r_iptr = (nk1-1) * siz_iz;
  MPI_Sendrecv(&x3d[s_iptr],1,DTypeZL,neighid[5],320,
               &x3d[r_iptr],1,DTypeZL,neighid[4],320,
               topocomm,&status);
  MPI_Sendrecv(&y3d[s_iptr],1,DTypeZL,neighid[5],320,
               &y3d[r_iptr],1,DTypeZL,neighid[4],320,
               topocomm,&status);
  MPI_Sendrecv(&z3d[s_iptr],1,DTypeZL,neighid[5],320,
               &z3d[r_iptr],1,DTypeZL,neighid[4],320,
               topocomm,&status);


  return 0;
}

int
grid_mesg_init(mympi_t *mympi, gd_t *gdcurv)
{
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;

  mympi->siz_sbuff_x1 = 3*nj*nk;
  mympi->siz_sbuff_x2 = 3*nj*nk;
  mympi->siz_sbuff_y1 = 3*ni*nk;
  mympi->siz_sbuff_y2 = 3*ni*nk;
  mympi->siz_sbuff_z1 = 3*ni*nj;
  mympi->siz_sbuff_z2 = 3*ni*nj;

  mympi->siz_rbuff_x1 = 3*nj*nk;
  mympi->siz_rbuff_x2 = 3*nj*nk;
  mympi->siz_rbuff_y1 = 3*ni*nk;
  mympi->siz_rbuff_y2 = 3*ni*nk;
  mympi->siz_rbuff_z1 = 3*ni*nj;
  mympi->siz_rbuff_z2 = 3*ni*nj;

  mympi->siz_sbuff = mympi->siz_sbuff_x1 + mympi->siz_sbuff_x2
                   + mympi->siz_sbuff_y1 + mympi->siz_sbuff_y2
                   + mympi->siz_sbuff_z1 + mympi->siz_sbuff_z2;

  mympi->siz_rbuff = mympi->siz_rbuff_x1 + mympi->siz_rbuff_x2
                   + mympi->siz_rbuff_y1 + mympi->siz_rbuff_y2
                   + mympi->siz_rbuff_z1 + mympi->siz_rbuff_z2;

  mympi->sbuff = (float *) malloc(mympi->siz_sbuff*sizeof(MPI_FLOAT));
  mympi->rbuff = (float *) malloc(mympi->siz_rbuff*sizeof(MPI_FLOAT));
  
  int tag[6] = {11, 12, 21, 22, 31, 32};

  mympi->s_reqs = (MPI_Request *) malloc(6*sizeof(MPI_Request));
  mympi->r_reqs = (MPI_Request *) malloc(6*sizeof(MPI_Request));
  // send
  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + mympi->siz_sbuff_x1;
  float *sbuff_y1 = sbuff_x2 + mympi->siz_sbuff_x2;
  float *sbuff_y2 = sbuff_y1 + mympi->siz_sbuff_y1;
  float *sbuff_z1 = sbuff_y2 + mympi->siz_sbuff_y2;
  float *sbuff_z2 = sbuff_z1 + mympi->siz_sbuff_z1;
  MPI_Send_init(sbuff_x1, mympi->siz_sbuff_x1, MPI_FLOAT, mympi->neighid[0], tag[0],
                mympi->topocomm, &(mympi->s_reqs[0]));
  MPI_Send_init(sbuff_x2, mympi->siz_sbuff_x2, MPI_FLOAT, mympi->neighid[1], tag[1],
                mympi->topocomm, &(mympi->s_reqs[1]));
  MPI_Send_init(sbuff_y1, mympi->siz_sbuff_y1, MPI_FLOAT, mympi->neighid[2], tag[2],
                mympi->topocomm, &(mympi->s_reqs[2]));
  MPI_Send_init(sbuff_y2, mympi->siz_sbuff_y2, MPI_FLOAT, mympi->neighid[3], tag[3],
                mympi->topocomm, &(mympi->s_reqs[3]));
  MPI_Send_init(sbuff_z1, mympi->siz_sbuff_z1, MPI_FLOAT, mympi->neighid[4], tag[4],
                mympi->topocomm, &(mympi->s_reqs[4]));
  MPI_Send_init(sbuff_z2, mympi->siz_sbuff_z2, MPI_FLOAT, mympi->neighid[5], tag[5],
                mympi->topocomm, &(mympi->s_reqs[5]));

  // recv
  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + mympi->siz_rbuff_x1;
  float *rbuff_y1 = rbuff_x2 + mympi->siz_rbuff_x2;
  float *rbuff_y2 = rbuff_y1 + mympi->siz_rbuff_y1;
  float *rbuff_z1 = rbuff_y2 + mympi->siz_rbuff_y2;
  float *rbuff_z2 = rbuff_z1 + mympi->siz_rbuff_z1;
  MPI_Recv_init(rbuff_x1, mympi->siz_rbuff_x1, MPI_FLOAT, mympi->neighid[0], tag[1],
                mympi->topocomm, &(mympi->r_reqs[0]));
  MPI_Recv_init(rbuff_x2, mympi->siz_rbuff_x2, MPI_FLOAT, mympi->neighid[1], tag[0],
                mympi->topocomm, &(mympi->r_reqs[1]));
  MPI_Recv_init(rbuff_y1, mympi->siz_rbuff_y1, MPI_FLOAT, mympi->neighid[2], tag[3],
                mympi->topocomm, &(mympi->r_reqs[2]));
  MPI_Recv_init(rbuff_y2, mympi->siz_rbuff_y2, MPI_FLOAT, mympi->neighid[3], tag[2],
                mympi->topocomm, &(mympi->r_reqs[3]));
  MPI_Recv_init(rbuff_z1, mympi->siz_rbuff_z1, MPI_FLOAT, mympi->neighid[4], tag[5],
                mympi->topocomm, &(mympi->r_reqs[4]));
  MPI_Recv_init(rbuff_z2, mympi->siz_rbuff_z2, MPI_FLOAT, mympi->neighid[5], tag[4],
                mympi->topocomm, &(mympi->r_reqs[5]));

  return 0;
}

int
grid_pack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x3d, float *y3d, float *z3d)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t iptr, iptr1;

  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + mympi->siz_sbuff_x1;
  float *sbuff_y1 = sbuff_x2 + mympi->siz_sbuff_x2;
  float *sbuff_y2 = sbuff_y1 + mympi->siz_sbuff_y1;
  float *sbuff_z1 = sbuff_y2 + mympi->siz_sbuff_y2;
  float *sbuff_z2 = sbuff_z1 + mympi->siz_sbuff_z1;

  // x1 
  for(int k=nk1; k<=nk2; k++) {
    for(int j=nj1; j<=nj2; j++)
    {
      iptr = k*siz_iz + j*siz_iy + ni1;
      iptr1 = (k-nk1)*nj + (j-nj1);
      sbuff_x1[iptr1]         = x3d[iptr];
      sbuff_x1[iptr1+nj*nk]   = y3d[iptr];
      sbuff_x1[iptr1+2*nj*nk] = z3d[iptr];
    }
  }

  // x2
  for(int k=nk1; k<=nk2; k++) {
    for(int j=nj1; j<=nj2; j++)
    {
      iptr = k*siz_iz + j*siz_iy + ni2;
      iptr1 = (k-nk1)*nj + (j-nj1);
      sbuff_x2[iptr1]         = x3d[iptr];
      sbuff_x2[iptr1+nj*nk]   = y3d[iptr];
      sbuff_x2[iptr1+2*nj*nk] = z3d[iptr];
    }
  }

  // y1 
  for(int k=nk1; k<=nk2; k++) {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = k*siz_iz + nj1*siz_iy + i;
      iptr1 = (k-nk1)*ni + (i-ni1);
      sbuff_y1[iptr1]         = x3d[iptr];
      sbuff_y1[iptr1+ni*nk]   = y3d[iptr];
      sbuff_y1[iptr1+2*ni*nk] = z3d[iptr];
    }
  }

  // y2 
  for(int k=nk1; k<=nk2; k++) {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = k*siz_iz + nj2*siz_iy + i;
      iptr1 = (k-nk1)*ni + (i-ni1);
      sbuff_y2[iptr1]         = x3d[iptr];
      sbuff_y2[iptr1+ni*nk]   = y3d[iptr];
      sbuff_y2[iptr1+2*ni*nk] = z3d[iptr];
    }
  }

  // z1 
  for(int j=nj1; j<=nj2; j++) {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = nk1*siz_iz + j*siz_iy + i;
      iptr1 = (j-nj1)*ni + (i-ni1);
      sbuff_z1[iptr1]         = x3d[iptr];
      sbuff_z1[iptr1+ni*nj]   = y3d[iptr];
      sbuff_z1[iptr1+2*ni*nj] = z3d[iptr];
    }
  }

  // z2 
  for(int j=nj1; j<=nj2; j++) {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = nk2*siz_iz + j*siz_iy + i;
      iptr1 = (j-nj1)*ni + (i-ni1);
      sbuff_z2[iptr1]         = x3d[iptr];
      sbuff_z2[iptr1+ni*nj]   = y3d[iptr];
      sbuff_z2[iptr1+2*ni*nj] = z3d[iptr];
    }
  }

  return 0;
}

int
grid_unpack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x3d, float *y3d, float *z3d)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int *neighid = mympi->neighid;
  size_t iptr, iptr1;

  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + mympi->siz_rbuff_x1;
  float *rbuff_y1 = rbuff_x2 + mympi->siz_rbuff_x2;
  float *rbuff_y2 = rbuff_y1 + mympi->siz_rbuff_y1;
  float *rbuff_z1 = rbuff_y2 + mympi->siz_rbuff_y2;
  float *rbuff_z2 = rbuff_z1 + mympi->siz_rbuff_z1;

  // x1 
  if(neighid[0] != MPI_PROC_NULL)
  {
    for(int k=nk1; k<=nk2; k++) {
      for(int j=nj1; j<=nj2; j++)
      {
        iptr = k*siz_iz + j*siz_iy + (ni1-1);
        iptr1 = (k-nk1)*nj + (j-nj1);
        x3d[iptr] = rbuff_x1[iptr1];
        y3d[iptr] = rbuff_x1[iptr1+nj*nk];
        z3d[iptr] = rbuff_x1[iptr1+2*nj*nk];
      }
    }
  }

  // x2 
  if(neighid[1] != MPI_PROC_NULL)
  {
    for(int k=nk1; k<=nk2; k++) {
      for(int j=nj1; j<=nj2; j++)
      {
        iptr = k*siz_iz + j*siz_iy + (ni2+1);
        iptr1 = (k-nk1)*nj + (j-nj1);
        x3d[iptr] = rbuff_x2[iptr1];
        y3d[iptr] = rbuff_x2[iptr1+nj*nk];
        z3d[iptr] = rbuff_x2[iptr1+2*nj*nk];
      }
    }
  }

  // y1 
  if(neighid[2] != MPI_PROC_NULL)
  {
    for(int k=nk1; k<=nk2; k++) {
      for(int i=ni1; i<=ni2; i++)
      {
        iptr = k*siz_iz + (nj1-1)*siz_iy + i;
        iptr1 = (k-nk1)*ni + (i-ni1);
        x3d[iptr] = rbuff_y1[iptr1];
        y3d[iptr] = rbuff_y1[iptr1+ni*nk];
        z3d[iptr] = rbuff_y1[iptr1+2*ni*nk];
      }
    }
  }

  // y2
  if(neighid[3] != MPI_PROC_NULL)
  {
    for(int k=nk1; k<=nk2; k++) {
      for(int i=ni1; i<=ni2; i++)
      {
        iptr = k*siz_iz + (nj2+1)*siz_iy + i;
        iptr1 = (k-nk1)*ni + (i-ni1);
        x3d[iptr] = rbuff_y2[iptr1];
        y3d[iptr] = rbuff_y2[iptr1+ni*nk];
        z3d[iptr] = rbuff_y2[iptr1+2*ni*nk];
      }
    }
  }

  // z1 
  if(neighid[4] != MPI_PROC_NULL)
  {
    for(int j=nj1; j<=nj2; j++) {
      for(int i=ni1; i<=ni2; i++)
      {
        iptr = (nk1-1)*siz_iz + j*siz_iy + i;
        iptr1 = (j-nj1)*ni + (i-ni1);
        x3d[iptr] = rbuff_z1[iptr1];
        y3d[iptr] = rbuff_z1[iptr1+ni*nj];
        z3d[iptr] = rbuff_z1[iptr1+2*ni*nj];
      }
    }
  }

  // z2 
  if(neighid[5] != MPI_PROC_NULL)
  {
    for(int j=nj1; j<=nj2; j++) {
      for(int i=ni1; i<=ni2; i++)
      {
        iptr = (nk2+1)*siz_iz + j*siz_iy + i;
        iptr1 = (j-nj1)*ni + (i-ni1);
        x3d[iptr] = rbuff_z2[iptr1];
        y3d[iptr] = rbuff_z2[iptr1+ni*nj];
        z3d[iptr] = rbuff_z2[iptr1+2*ni*nj];
      }
    }
  }

  return 0;
}
