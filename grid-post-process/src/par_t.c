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

  par->num_of_grid = 0;
  if (item = cJSON_GetObjectItem(root, "input_grids_info")) 
  {
    par->num_of_grid = cJSON_GetArraySize(item);
    par->num_of_points = (int *)malloc(par->num_of_grid*sizeof(int)*CONST_NDIM);
    par->num_of_procs_in = (int *)malloc(par->num_of_grid*sizeof(int)*CONST_NDIM);
    par->flag_stretch = (int *)malloc(par->num_of_grid*sizeof(int));
    par->import_dir = (char **)malloc(par->num_of_grid*sizeof(char*));
    par->stretch_file = (char **)malloc(par->num_of_grid*sizeof(char*));
    // each input grid info
    for (int i=0; i < par->num_of_grid; i++)
    {
      par->import_dir[i] = (char *)malloc(PAR_MAX_STRLEN*sizeof(char));
      subitem = cJSON_GetArrayItem(item, i);
      if (thirditem = cJSON_GetObjectItem(subitem, "grid_import_dir"))
      {
        sprintf(par->import_dir[i],"%s",thirditem->valuestring);
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "number_of_grid_points"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->num_of_points[i*CONST_NDIM+j] = cJSON_GetArrayItem(thirditem, j)->valueint;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "number_of_mpiprocs_in"))
      {
        for (int j = 0; j < CONST_NDIM; j++) {
          par->num_of_procs_in[i*CONST_NDIM+j] = cJSON_GetArrayItem(thirditem, j)->valueint;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "flag_stretch"))
      {
        par->flag_stretch[i] = thirditem->valueint;
      }
      if(par->flag_stretch[i] == 1)
      {
        par->stretch_file[i] = (char *)malloc(PAR_MAX_STRLEN*sizeof(char));
        if (thirditem = cJSON_GetObjectItem(subitem, "stretch_file"))
        {
          sprintf(par->stretch_file[i],"%s",thirditem->valuestring);
        }
      }
    }
  }

  par->num_of_procs_out = (int *)malloc(sizeof(int)*CONST_NDIM);
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_out")) 
  {
    for (int i=0; i < CONST_NDIM; i++)
    {
      par->num_of_procs_out[i] = cJSON_GetArrayItem(item, i)->valueint;
    }
  }

  par->stretch_idire = 0;
  if (item = cJSON_GetObjectItem(root, "stretch_direction")) 
  {
    sprintf(par->stretch_dire, "%s", item->valuestring);
    if(strcmp(par->stretch_dire, "x")==0)
    {
      par->stretch_idire = X_DIRE;
    } else if (strcmp(par->stretch_dire, "y")==0) {
      par->stretch_idire = Y_DIRE;
    } else if (strcmp(par->stretch_dire, "z")==0) {
      par->stretch_idire = Z_DIRE;
    } else {
      fprintf(stdout,"ERROR: must give stretch direction\n");
      exit(-1);
    }
  }
  par->merge_idire = 0;
  if (par->num_of_grid > 1)
  {
    if (item = cJSON_GetObjectItem(root, "merge_direction")) 
    {
      sprintf(par->merge_dire, "%s", item->valuestring);
      if(strcmp(par->merge_dire, "x")==0)
      {
        par->merge_idire = X_DIRE;
      } else if (strcmp(par->merge_dire, "y")==0) {
        par->merge_idire = Y_DIRE;
      } else if (strcmp(par->merge_dire, "z")==0) {
        par->merge_idire = Z_DIRE;
      }
    }
    if(par->merge_idire == 0)
    {
      fprintf(stdout,"ERROR: gird number > 1, must give merge direction\n");
      exit(1);
    }
  }

  // default pml layers
  par->flag_pml = 0;
  par->num_of_pml_x1 = 0;
  par->num_of_pml_x2 = 0;
  par->num_of_pml_y1 = 0;
  par->num_of_pml_y2 = 0;
  par->num_of_pml_z1 = 0;
  par->num_of_pml_z2 = 0;

  // flag_pml = 1; double pml layer for computation balance
  // eg. nx = 100, pml_x1 = 10, pml_x2 = 10, nproc_x=3
  // gird points is (30 40 30)
  // not (33,33,34)
  if (item = cJSON_GetObjectItem(root, "flag_pml")) {
    par->flag_pml = item->valueint;
  }
  if(par->flag_pml==1)
  {
    if (item = cJSON_GetObjectItem(root, "pml_layers")) {
      if (subitem = cJSON_GetObjectItem(item, "number_of_pml_x1")) {
        par->num_of_pml_x1 = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(item, "number_of_pml_x2")) {
        par->num_of_pml_x2 = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(item, "number_of_pml_y1")) {
        par->num_of_pml_y1 = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(item, "number_of_pml_y2")) {
        par->num_of_pml_y2 = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(item, "number_of_pml_z1")) {
        par->num_of_pml_z1 = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(item, "number_of_pml_z2")) {
        par->num_of_pml_z2 = subitem->valueint;
      }
    }
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

  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
    sprintf(par->export_dir, "%s", item->valuestring);
  }

  return ierr;
}


int
par_print(par_t *par)
{    
  int ierr = 0;
  
  for(int i=0; i<par->num_of_grid; i++)
  {
    fprintf(stdout,"input grid index number is %d\n",i+1);

    fprintf(stdout,"number of total gird points is [%d,%d,%d]\n",
    par->num_of_points[0+i*CONST_NDIM],par->num_of_points[1+i*CONST_NDIM],
    par->num_of_points[2+i*CONST_NDIM]);

    fprintf(stdout,"number of mpi procs input is [%d,%d,%d]\n",
    par->num_of_procs_in[0+i*CONST_NDIM],par->num_of_procs_in[1+i*CONST_NDIM],
    par->num_of_procs_in[2+i*CONST_NDIM]);

    fprintf(stdout,"input grid dir is \n %s\n", par->import_dir[i]);
    if(par->flag_stretch[i] == 1) {
      fprintf(stdout, "------- strech direction is %s\n",par->stretch_dire);
      fprintf(stdout, "------- strech file is %s\n",par->stretch_file[i]);
    }
  }

  if(par->num_of_grid>1)
  {
    fprintf(stdout,"merge direction is %s\n",par->merge_dire);
  }

  fprintf(stdout,"number of mpi procs output is [%d,%d,%d]\n",
  par->num_of_procs_out[0], par->num_of_procs_out[1],
  par->num_of_procs_out[2]);

  fprintf(stdout,"number of pml layers x1 is %d\n",par->num_of_pml_x1);
  fprintf(stdout,"number of pml layers x2 is %d\n",par->num_of_pml_x2);
  fprintf(stdout,"number of pml layers y1 is %d\n",par->num_of_pml_y1);
  fprintf(stdout,"number of pml layers y2 is %d\n",par->num_of_pml_y2);
  fprintf(stdout,"number of pml layers z1 is %d\n",par->num_of_pml_z1);
  fprintf(stdout,"number of pml layers z2 is %d\n",par->num_of_pml_z2);

  fprintf(stdout, "-------------------------------------------------------\n");
  if (par->flag_sample == 1) {
    fprintf(stdout,"------- sample grid xi direction factor is %d------- \n",par->sample_factor_xi);
    fprintf(stdout,"------- sample grid et direction factor is %d------- \n",par->sample_factor_et);
    fprintf(stdout,"------- sample grid zt direction factor is %d------- \n",par->sample_factor_zt);
  }

  fprintf(stdout,"export grid dir is \n %s\n",par->export_dir);
  fflush(stdout);

  return ierr;
}
