#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "par_t.h"

/*
 * for MPI, master read, broadcast to all procs
 */
int
par_mpi_get(char *par_fname, int myid, MPI_Comm comm, par_t *par)
{
  char *str;

  if (myid==0)
  {
    FILE *fp=fopen(par_fname,"r");
    if (!fp) {
      fprintf(stderr,"Error: can't open par file: %s\n", par_fname);
      MPI_Abort(MPI_COMM_WORLD,9);
    }
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);

    // bcast len to all nodes
    MPI_Bcast(&len, 1, MPI_LONG, 0, comm);

    fseek(fp, 0, SEEK_SET);
    str = (char*)malloc(len+1);
    fread(str, 1, len, fp);
    fclose(fp);

    // bcast str
    MPI_Bcast(str, len+1, MPI_CHAR, 0, comm);
  }
  else
  {
    long len;
    // get len 
    MPI_Bcast(&len, 1, MPI_LONG, 0, comm);

    str = (char*)malloc(len+1);
    // get str
    MPI_Bcast(str, len+1, MPI_CHAR, 0, comm);
  }
    
  par_read_from_str(str, par);

  free(str);

  return 0;
}

/*
 * funcs to get par from alread read in str
 */
int 
par_read_from_str(const char *str, par_t *par)
{
  int ierr = 0;

  // convert str to json
  cJSON *root = cJSON_Parse(str);
  if (NULL == root) {
    printf("Error at parsing json!\n");
    exit(1);
  }

  cJSON *item;
  cJSON *subitem, *thirditem;

  // default mpi threads
  par->number_of_mpiprocs_x = 1;
  par->number_of_mpiprocs_y = 1;
  par->number_of_mpiprocs_z = 1;

  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_x")) {
    par->number_of_mpiprocs_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_y")) {
    par->number_of_mpiprocs_y = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_z")) {
    par->number_of_mpiprocs_z = item->valueint;
  }

  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_x")) {
    par->number_of_grid_points_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_y")) {
    par->number_of_grid_points_y = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_z")) {
    par->number_of_grid_points_z = item->valueint;
  }

  // default not check
  par->grid_check  = 0;
  par->check_orth  = 0;
  par->check_jac   = 0;
  par->check_step_xi  = 0;
  par->check_step_et  = 0;
  par->check_step_zt  = 0;
  par->check_smooth_xi = 0;
  par->check_smooth_et = 0;
  par->check_smooth_zt = 0;
  if (item = cJSON_GetObjectItem(root, "check_orth")) {
    par->check_orth = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_jac")) {
    par->check_jac = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_xi")) {
    par->check_step_xi = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_et")) {
    par->check_step_et = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_zt")) {
    par->check_step_zt = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_xi")) {
    par->check_smooth_xi = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_et")) {
    par->check_smooth_et = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_zt")) {
    par->check_smooth_zt = item->valueint;
  }
  int check = par->check_orth + par->check_jac
            + par->check_step_xi + par->check_step_et
            + par->check_step_zt + par->check_smooth_xi
            + par->check_smooth_et + par->check_smooth_zt;
  if(check != 0)
  {
    par->grid_check = 1;
  }

  if (item = cJSON_GetObjectItem(root, "geometry_input_file")) {
    sprintf(par->geometry_input_file, "%s", item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
    sprintf(par->output_dir, "%s", item->valuestring);
  }
 
  if (item = cJSON_GetObjectItem(root, "grid_method")) {
    if (subitem = cJSON_GetObjectItem(item, "linear_tfi")) {
      par->method_itype = LINEAR_TFI;
    }
    if (subitem = cJSON_GetObjectItem(item, "elli_diri")) {
      par->method_itype = ELLI_DIRI;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        for (int i=0; i<6; i++)
        {
          par->coef[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "iter_err")) {
        par->iter_err = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "max_iter")) {
        par->max_iter = thirditem->valueint;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "elli_higen")) {
      par->method_itype = ELLI_HIGEN;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        for (int i=0; i<6; i++)
        {
          par->coef[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "iter_err")) {
        par->iter_err = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "max_iter")) {
        par->max_iter = thirditem->valueint;
      }
    }
  }

  return ierr;
}


int
par_print(par_t *par)
{    
  int ierr = 0;

  fprintf(stdout,"number of total gird points x is %d\n",par->number_of_grid_points_x);
  fprintf(stdout,"number of total gird points y is %d\n",par->number_of_grid_points_y);
  fprintf(stdout,"number of total gird points z is %d\n",par->number_of_grid_points_z);

  fprintf(stdout,"number of mpi procs x is %d\n",par->number_of_mpiprocs_x);
  fprintf(stdout,"number of mpi procs y is %d\n",par->number_of_mpiprocs_y);
  fprintf(stdout,"number of mpi procs z is %d\n",par->number_of_mpiprocs_z);

  fprintf(stdout,"input geometry file is \n %s\n",par->geometry_input_file);
  fprintf(stdout,"export grid dir is \n %s\n",par->output_dir);
  fprintf(stdout, "-------------------------------------------------------\n");
  if (par->grid_check == 1) {
    fprintf(stdout, "------- grid quality check-------\n");
  }
  if (par->check_orth == 1) {
    fprintf(stdout, "------- check grid orthogonality-------\n");
  }
  if (par->check_jac == 1) {
    fprintf(stdout, "------- check grid jacobi-------\n");
  }
  if (par->check_step_xi == 1) {
    fprintf(stdout, "------- check grid step xi direction-------\n");
  }
  if (par->check_step_et == 1) {
    fprintf(stdout, "------- check grid step et direction-------\n");
  }
  if (par->check_step_zt == 1) {
    fprintf(stdout, "------- check grid step zt direction-------\n");
  }
  if (par->check_smooth_xi == 1) {
    fprintf(stdout, "------- check grid smooth xi direction-------\n");
  }
  if (par->check_smooth_et == 1) {
    fprintf(stdout, "------- check grid smooth et direction-------\n");
  }
  if (par->check_smooth_zt == 1) {
    fprintf(stdout, "------- check grid smooth zt direction-------\n");
  }

  fprintf(stdout, "------- grid generate method-------\n");
  if(par->method_itype == LINEAR_TFI) {
    fprintf(stdout, "grid generate method is linear_tfi\n");
  }
  if(par->method_itype == ELLI_DIRI) {
    fprintf(stdout, "grid generate method is elliptic_dirichlet\n");
    fprintf(stdout, "elli_diri x1 coef is %f\n", par->coef[0]);
    fprintf(stdout, "elli_diri x2 coef is %f\n", par->coef[1]);
    fprintf(stdout, "elli_diri y1 coef is %f\n", par->coef[2]);
    fprintf(stdout, "elli_diri y2 coef is %f\n", par->coef[3]);
    fprintf(stdout, "elli_diri z1 coef is %f\n", par->coef[4]);
    fprintf(stdout, "elli_diri z2 coef is %f\n", par->coef[5]);
    fprintf(stdout, "max_iteration is %d\n", par->max_iter);
    fprintf(stdout, "iter_error is %f\n", par->iter_err);
  }
  if(par->method_itype == ELLI_HIGEN) {
    fprintf(stdout, "grid generate method is elliptic_hilgenstock\n");
    fprintf(stdout, "elli_higen x1 coef is %f\n", par->coef[0]);
    fprintf(stdout, "elli_higen x2 coef is %f\n", par->coef[1]);
    fprintf(stdout, "elli_higen y1 coef is %f\n", par->coef[2]);
    fprintf(stdout, "elli_higen y2 coef is %f\n", par->coef[3]);
    fprintf(stdout, "elli_higen z1 coef is %f\n", par->coef[4]);
    fprintf(stdout, "elli_higen z2 coef is %f\n", par->coef[5]);
    fprintf(stdout, "max_iteration is %d\n", par->max_iter);
    fprintf(stdout, "iter_error is %f\n", par->iter_err);
  }

  return ierr;
}
