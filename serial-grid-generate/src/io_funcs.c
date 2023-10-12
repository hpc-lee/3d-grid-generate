
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
  int nprocx_out = par->number_of_mpiprocs_x_out;
  int nprocy_out = par->number_of_mpiprocs_y_out;
  int nprocz_out = par->number_of_mpiprocs_z_out;

  int total_nx = gdcurv->nx;
  int total_ny = gdcurv->ny;
  int total_nz = gdcurv->nz;

  int gni1, gnj1, gnk1;
  int gni, gnj, gnk;
  int ni, nj, nk;
  size_t iptr, iptr1;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;
  float *coord_x = (float *) malloc(sizeof(float)*total_nx*total_ny*total_nz);
  float *coord_y = (float *) malloc(sizeof(float)*total_nx*total_ny*total_nz);
  float *coord_z = (float *) malloc(sizeof(float)*total_nx*total_ny*total_nz);

  // construct file name
  char fname_coords[CONST_MAX_STRLEN];
  char ou_file[CONST_MAX_STRLEN];
  
  // read in nc
  int ierr;
  int ncid;
  int xid, yid, zid;
  int dimid[3];
  int g_start[3];
  int count[3];

  // export coord nc file
  for(int kk=0; kk<nprocz_out; kk++) {
    for(int jj=0; jj<nprocy_out; jj++) {
      for(int ii=0; ii<nprocx_out; ii++)
      {
        sprintf(fname_coords,"px%d_py%d_pz%d",ii,jj,kk);
        sprintf(ou_file,"%s/coord_%s.nc",par->grid_export_dir,fname_coords);

        gd_info_set(par, ii, jj, kk, g_start, count);
        gni1 = g_start[0];
        gnj1 = g_start[1];
        gnk1 = g_start[2];
        ni = count[0];
        nj = count[1];
        nk = count[2];
        
        for(int k=0; k<nk; k++) {
          for(int j=0; j<nj; j++) {
            for(int i=0; i<ni; i++)
            {
              gni = gni1 + i;
              gnj = gnj1 + j;
              gnk = gnk1 + k;
              iptr = gnk*total_nx*total_ny + gnj*total_nx + gni;

              iptr1 = k*ni*nj + j*ni + i;

              coord_x[iptr1] = x3d[iptr];
              coord_y[iptr1] = y3d[iptr];
              coord_z[iptr1] = z3d[iptr];
            }
          }
        }

        ierr = nc_create(ou_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid); handle_nc_err(ierr);

        // define dimension
        ierr = nc_def_dim(ncid, "i", ni, &dimid[2]); handle_nc_err(ierr);
        ierr = nc_def_dim(ncid, "j", nj, &dimid[1]); handle_nc_err(ierr);
        ierr = nc_def_dim(ncid, "k", nk, &dimid[0]); handle_nc_err(ierr);

        // define vars
        ierr = nc_def_var(ncid, "x", NC_FLOAT, 3, dimid, &xid); handle_nc_err(ierr);
        ierr = nc_def_var(ncid, "y", NC_FLOAT, 3, dimid, &yid); handle_nc_err(ierr);
        ierr = nc_def_var(ncid, "z", NC_FLOAT, 3, dimid, &zid); handle_nc_err(ierr);

        // attribute: index in output snapshot, index w ghost in thread
        nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                         NC_INT,3,g_start);

        nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                         NC_INT,3,count);

        // end def
        ierr = nc_enddef(ncid);
        handle_nc_err(ierr);

        // add vars
        ierr = nc_put_var_float(ncid, xid, coord_x);  handle_nc_err(ierr);
        ierr = nc_put_var_float(ncid, yid, coord_y);  handle_nc_err(ierr);
        ierr = nc_put_var_float(ncid, zid, coord_z);  handle_nc_err(ierr);
        
        // close file
        ierr = nc_close(ncid);
        handle_nc_err(ierr);
      }
    }
  }

  free(coord_x);
  free(coord_y);
  free(coord_z);

  return 0;
}

int
quality_export(io_quality_t *io_quality, par_t *par, char *var_name)
{
  int nprocx_out = par->number_of_mpiprocs_x_out;
  int nprocy_out = par->number_of_mpiprocs_y_out;
  int nprocz_out = par->number_of_mpiprocs_z_out;

  int total_nx = io_quality->nx;
  int total_ny = io_quality->ny;
  int total_nz = io_quality->nz;

  int gni1, gnj1, gnk1;
  int gni, gnj, gnk;
  int ni, nj, nk;
  size_t iptr, iptr1;

  float *var = io_quality->var;
  float *var_out = (float *) malloc(sizeof(float)*total_nx*total_ny*total_nz);

  // construct file name
  char fname_coords[CONST_MAX_STRLEN];
  char ou_file[CONST_MAX_STRLEN];
  
  // read in nc
  int ierr;
  int ncid;
  int varid;
  int dimid[3];
  int g_start[3];
  int count[3];

  // export coord nc file
  for(int kk=0; kk<nprocz_out; kk++) {
    for(int jj=0; jj<nprocy_out; jj++) {
      for(int ii=0; ii<nprocx_out; ii++)
      {
        sprintf(fname_coords,"px%d_py%d_pz%d",ii,jj,kk);
        sprintf(ou_file,"%s/%s_%s.nc",par->grid_export_dir,var_name,fname_coords);

        gd_info_set(par, ii, jj, kk, g_start, count);
        gni1 = g_start[0];
        gnj1 = g_start[1];
        gnk1 = g_start[2];
        ni = count[0];
        nj = count[1];
        nk = count[2];
        
        for(int k=0; k<nk; k++) {
          for(int j=0; j<nj; j++) {
            for(int i=0; i<ni; i++)
            {
              gni = gni1 + i;
              gnj = gnj1 + j;
              gnk = gnk1 + k;
              iptr = gnk*total_nx*total_ny + gnj*total_nx + gni;

              iptr1 = k*ni*nj + j*ni + i;

              var_out[iptr1] = var[iptr];
            }
          }
        }

        ierr = nc_create(ou_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid); handle_nc_err(ierr);

        // define dimension
        ierr = nc_def_dim(ncid, "i", ni, &dimid[2]); handle_nc_err(ierr);
        ierr = nc_def_dim(ncid, "j", nj, &dimid[1]); handle_nc_err(ierr);
        ierr = nc_def_dim(ncid, "k", nk, &dimid[0]); handle_nc_err(ierr);

        // define vars
        ierr = nc_def_var(ncid, var_name, NC_FLOAT, 3, dimid, &varid); handle_nc_err(ierr);

        // attribute: index in output snapshot, index w ghost in thread
        nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                         NC_INT,3,g_start);

        nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                         NC_INT,3,count);

        // end def
        ierr = nc_enddef(ncid);
        handle_nc_err(ierr);

        // add vars
        ierr = nc_put_var_float(ncid, varid, var_out);  handle_nc_err(ierr);
        
        // close file
        ierr = nc_close(ncid);
        handle_nc_err(ierr);
      }
    }
  }

  free(var_out);

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
