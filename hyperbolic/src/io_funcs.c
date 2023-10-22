
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
gd_curv_coord_export(gd_t *gdcurv, par_t *par)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;

  size_t iptr, iptr1;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  // construct file name
  char fname_coords[CONST_MAX_STRLEN];
  char ou_file[CONST_MAX_STRLEN];
  
  // read in nc
  int ierr;
  int ncid;
  int xid, yid, zid;
  int dimid[3];

  // export coord nc file
  sprintf(fname_coords,"px%d_py%d_pz%d",0,0,0);
  sprintf(ou_file,"%s/coord_%s.nc",par->grid_export_dir,fname_coords);

  ierr = nc_create(ou_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid); handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]); handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]); handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]); handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, "x", NC_FLOAT, 3, dimid, &xid); handle_nc_err(ierr);
  ierr = nc_def_var(ncid, "y", NC_FLOAT, 3, dimid, &yid); handle_nc_err(ierr);
  ierr = nc_def_var(ncid, "z", NC_FLOAT, 3, dimid, &zid); handle_nc_err(ierr);

  // attribute: index in output snapshot, index w ghost in thread
  int g_start[3] = {0, 0, 0};
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,3,g_start);

  int count[3] = {nx, ny, nz};
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,3,count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  ierr = nc_put_var_float(ncid, xid, x3d);  handle_nc_err(ierr);
  ierr = nc_put_var_float(ncid, yid, y3d);  handle_nc_err(ierr);
  ierr = nc_put_var_float(ncid, zid, z3d);  handle_nc_err(ierr);
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return 0;
}

int
quality_export(io_quality_t *io_quality, par_t *par, char *var_name)
{
  int nx = io_quality->nx;
  int ny = io_quality->ny;
  int nz = io_quality->nz;

  size_t iptr, iptr1;

  float *var = io_quality->var;

  // construct file name
  char fname_coords[CONST_MAX_STRLEN];
  char ou_file[CONST_MAX_STRLEN];
  
  // read in nc
  int ierr;
  int ncid;
  int varid;
  int dimid[3];

  // export coord nc file
  sprintf(fname_coords,"px%d_py%d_pz%d",0,0,0);
  sprintf(ou_file,"%s/%s_%s.nc",par->grid_export_dir,var_name,fname_coords);
  
  ierr = nc_create(ou_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid); handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]); handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]); handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]); handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, var_name, NC_FLOAT, 3, dimid, &varid); handle_nc_err(ierr);

  // attribute: index in output snapshot, index w ghost in thread
  int g_start[3] = {0, 0, 0};
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,3,g_start);

  int count[3] = {nx, ny, nz};
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,3,count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  ierr = nc_put_var_float(ncid, varid, var);  handle_nc_err(ierr);
  
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
