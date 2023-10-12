#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "par_t.h"

/*
 * read from file
 */

int
par_read_from_file(char *par_fname,  par_t *par, int verbose)
{
  //
  // read whole file into str
  //
  FILE *fp = fopen(par_fname,"r");
  if (!fp) {
    fprintf(stderr,"Error: can't open par file: %s\n", par_fname);
    exit(1);
  }

  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);

  fseek(fp, 0, SEEK_SET);
  char *str = (char*)malloc(len+1);
  fread(str, 1, len, fp);
  fclose(fp);

  // read from str
  par_read_from_str(str, par);

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

  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_x")) {
    par->number_of_grid_points_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_y")) {
    par->number_of_grid_points_y = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_z")) {
    par->number_of_grid_points_z = item->valueint;
  }

  // default mpi threads
  par->number_of_mpiprocs_x_in = 1;
  par->number_of_mpiprocs_y_in = 1;
  par->number_of_mpiprocs_z_in = 1;

  par->number_of_mpiprocs_x_out = 1;
  par->number_of_mpiprocs_y_out = 1;
  par->number_of_mpiprocs_z_out = 1;

  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_x_in")) {
    par->number_of_mpiprocs_x_in = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_y_in")) {
    par->number_of_mpiprocs_y_in = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_z_in")) {
    par->number_of_mpiprocs_z_in = item->valueint;
  }

  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_x_out")) {
    par->number_of_mpiprocs_x_out = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_y_out")) {
    par->number_of_mpiprocs_y_out = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_z_out")) {
    par->number_of_mpiprocs_z_out = item->valueint;
  }

  // default pml layers
  par->number_of_pml_x1 = 0;
  par->number_of_pml_x2 = 0;
  par->number_of_pml_y1 = 0;
  par->number_of_pml_y2 = 0;
  par->number_of_pml_z1 = 0;
  par->number_of_pml_z2 = 0;
  if (item = cJSON_GetObjectItem(root, "pml_layers")) {
    if (subitem = cJSON_GetObjectItem(item, "number_of_pml_x1")) {
      par->number_of_pml_x1 = subitem->valueint;
    }
    if (subitem = cJSON_GetObjectItem(item, "number_of_pml_x2")) {
      par->number_of_pml_x2 = subitem->valueint;
    }
    if (subitem = cJSON_GetObjectItem(item, "number_of_pml_y1")) {
      par->number_of_pml_y1 = subitem->valueint;
    }
    if (subitem = cJSON_GetObjectItem(item, "number_of_pml_y2")) {
      par->number_of_pml_y2 = subitem->valueint;
    }
    if (subitem = cJSON_GetObjectItem(item, "number_of_pml_z1")) {
      par->number_of_pml_z1 = subitem->valueint;
    }
    if (subitem = cJSON_GetObjectItem(item, "number_of_pml_z2")) {
      par->number_of_pml_z2 = subitem->valueint;
    }
  }
  // default not strech
  par->flag_strech_xi = 0;
  par->flag_strech_et = 0;
  par->flag_strech_zt = 0;
  if (item = cJSON_GetObjectItem(root, "flag_strech_xi")) {
    par->flag_strech_xi = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_xi_coef")) {
    par->strech_xi_coef = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "flag_strech_et")) {
    par->flag_strech_et = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_et_coef")) {
    par->strech_et_coef = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "flag_strech_zt")) {
    par->flag_strech_zt = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_zt_coef")) {
    par->strech_zt_coef = item->valuedouble;
  }

  // default not intep
  par->sample_factor_xi = 1;
  par->sample_factor_et = 1;
  par->sample_factor_zt = 1;
  if (item = cJSON_GetObjectItem(root, "flag_sample")) {
    par->flag_sample = item->valueint;
  }
  if (par->flag_sample == 1) 
  {
    if (item = cJSON_GetObjectItem(root, "sample_factor_xi")) {
      par->sample_factor_xi = item->valueint;
      if((par->sample_factor_xi-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_xi must >= 1\n");
        exit(1);
      }
    }
    if (item = cJSON_GetObjectItem(root, "sample_factor_et")) {
      par->sample_factor_et = item->valueint;
      if((par->sample_factor_et-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_et must >= 1\n");
        exit(1);
      }
    }
    if (item = cJSON_GetObjectItem(root, "sample_factor_zt")) {
      par->sample_factor_zt = item->valueint;
      if((par->sample_factor_zt-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_zt must >= 1\n");
        exit(1);
      }
    }
  }

  if (item = cJSON_GetObjectItem(root, "grid_import_dir")) {
    sprintf(par->import_dir, "%s", item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
    sprintf(par->export_dir, "%s", item->valuestring);
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

  fprintf(stdout,"number of mpi procs x import is %d\n",par->number_of_mpiprocs_x_in);
  fprintf(stdout,"number of mpi procs y import is %d\n",par->number_of_mpiprocs_y_in);
  fprintf(stdout,"number of mpi procs z import is %d\n",par->number_of_mpiprocs_z_in);

  fprintf(stdout,"number of mpi procs x export is %d\n",par->number_of_mpiprocs_x_out);
  fprintf(stdout,"number of mpi procs y export is %d\n",par->number_of_mpiprocs_y_out);
  fprintf(stdout,"number of mpi procs z export is %d\n",par->number_of_mpiprocs_z_out);

  fprintf(stdout,"number of pml layers x1 is %d\n",par->number_of_pml_x1);
  fprintf(stdout,"number of pml layers x2 is %d\n",par->number_of_pml_x2);
  fprintf(stdout,"number of pml layers y1 is %d\n",par->number_of_pml_y1);
  fprintf(stdout,"number of pml layers y2 is %d\n",par->number_of_pml_y2);
  fprintf(stdout,"number of pml layers z1 is %d\n",par->number_of_pml_z1);
  fprintf(stdout,"number of pml layers z2 is %d\n",par->number_of_pml_z2);

  fprintf(stdout, "-------------------------------------------------------\n");
  if (par->flag_sample == 1) {
    fprintf(stdout,"------- sample grid xi direction factor is %d------- \n",par->sample_factor_xi);
    fprintf(stdout,"------- sample grid et direction factor is %d------- \n",par->sample_factor_et);
    fprintf(stdout,"------- sample grid zt direction factor is %d------- \n",par->sample_factor_zt);
  }
  if(par->flag_strech_xi == 1) {
    fprintf(stdout, "------- strech xi and strech coef is %f-------\n",par->strech_xi_coef);
  }
  if(par->flag_strech_et == 1) {
    fprintf(stdout, "------- strech et and strech coef is %f-------\n",par->strech_et_coef);
  }
  if(par->flag_strech_zt == 1) {
    fprintf(stdout, "------- strech zt and strech coef is %f-------\n",par->strech_zt_coef);
  }

  fprintf(stdout,"input grid dir is \n %s\n", par->import_dir);
  fprintf(stdout,"export grid dir is \n %s\n",par->export_dir);
  fflush(stdout);

  return ierr;
}
