
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h" 

#include "lib_mem.h"
#include "constants.h"
#include "io_funcs.h"

int
init_io_quality(io_quality_t *io_quality, gd_t *gdcurv)
{
  io_quality->nx = gdcurv->nx;
  io_quality->ny = gdcurv->ny;
  io_quality->nz = gdcurv->nz;
  
  // malloc quality space 
  io_quality->var = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp, 0.0, "quality_init");
  if (io_quality->var == NULL) {
      fprintf(stderr,"Error: failed to alloc quality vars\n");
      fflush(stderr);
  }
   
  return 0;
}

int
gd_curv_coord_export(gd_t *gdcurv, char *output_dir)
{
  size_t *restrict c3d_pos   = gdcurv->cmp_pos;
  char  **restrict c3d_name  = gdcurv->cmp_name;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[gdcurv->ncmp];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    ierr = nc_def_var(ncid, gdcurv->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
    handle_nc_err(ierr);
  }

  int l_count[] = { nx, ny, nz };
  nc_put_att_int(ncid,NC_GLOBAL,"number_of_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    float *ptr = gdcurv->v4d + gdcurv->cmp_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return 0;
}

int
quality_export(io_quality_t *io_quality, char *output_dir, char *var_name)
{
  int  nx = io_quality->nx;
  int  ny = io_quality->ny;
  int  nz = io_quality->nz;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/%s.nc", output_dir,var_name);
  
  // read in nc
  int ncid;
  int varid;
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, var_name, NC_FLOAT, CONST_NDIM, dimid, &varid);
  handle_nc_err(ierr);

  int l_count[] = { nx, ny, nz };
  nc_put_att_int(ncid,NC_GLOBAL,"number_of_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  float *ptr = io_quality->var;
  ierr = nc_put_var_float(ncid,varid,ptr);
  handle_nc_err(ierr);
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return 0;
}

int
io_get_nextline(FILE *fp, char *str, int length)
{
  int ierr = 0;

  do
  {
    if (fgets(str, length, fp) == NULL)
    {
      ierr = 1;
      return ierr;
    }
  } while (str[0] == '#' || str[0] == '\n');

  // remove newline char
  int len = strlen(str);

  if (len > 0 && str[len-1] == '\n') {
    str[len-1] = '\0';
  }

  // for debug:
  //fprintf(stdout," --return: %s\n", str);

  return ierr;
}
