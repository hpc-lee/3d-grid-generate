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

  // default not check
  par->grid_check = 0;
  par->check_orth  = 0;
  par->check_jac   = 0;
  par->check_ratio = 0;
  par->check_step_x  = 0;
  par->check_step_z  = 0;
  par->check_smooth_x = 0;
  par->check_smooth_z = 0;
  if (item = cJSON_GetObjectItem(root, "check_orth")) {
    par->check_orth = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_jac")) {
    par->check_jac = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_ratio")) {
    par->check_ratio = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_x")) {
    par->check_step_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_z")) {
    par->check_step_z = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_x")) {
    par->check_smooth_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_z")) {
    par->check_smooth_z = item->valueint;
  }
  int check = par->check_orth + par->check_jac /
            + par->check_ratio + par->check_step_x /
            + par->check_step_z +  par->check_smooth_x /
            + par->check_smooth_z; 
  if(check != 0)
  {
    par->grid_check = 1;
  }

  // default not strech
  par->flag_strech_x = 0;
  par->flag_strech_z = 0;
  if (item = cJSON_GetObjectItem(root, "flag_strech_x")) {
    par->flag_strech_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_x_coef")) {
    par->strech_x_coef = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "flag_strech_z")) {
    par->flag_strech_z = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_z_coef")) {
    par->strech_z_coef = item->valuedouble;
  }

  // default not intep
  par->sample_factor_x = 1.0;
  par->sample_factor_z = 1.0;
  if (item = cJSON_GetObjectItem(root, "flag_sample_x")) {
    par->flag_sample_x = item->valueint;
  }
  if (par->flag_sample_x == 1) 
  {
    if (item = cJSON_GetObjectItem(root, "sample_factor_x")) {
      par->sample_factor_x = item->valuedouble;
      if((par->sample_factor_x-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_x must >= 1\n");
        exit(1);
      }
    }
  }
  if (item = cJSON_GetObjectItem(root, "flag_sample_z")) {
    par->flag_sample_z = item->valueint;
  }
  if (par->flag_sample_z == 1) 
  {
    if (item = cJSON_GetObjectItem(root, "sample_factor_z")) {
      par->sample_factor_z = item->valuedouble;
      if((par->sample_factor_z-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_z must >= 1\n");
        exit(1);
      }
    }
  }

  if (item = cJSON_GetObjectItem(root, "geometry_input_file")) {
    sprintf(par->geometry_input_file, "%s", item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
    sprintf(par->grid_export_dir, "%s", item->valuestring);
  }

  if (item = cJSON_GetObjectItem(root, "grid_method")) {
    if (subitem = cJSON_GetObjectItem(item, "linear_TFI")) {
      par->method_itype = TFI;
    }
    if (subitem = cJSON_GetObjectItem(item, "hermite")) {
      par->method_itype = HERMITE;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "elli_diri")) {
      par->method_itype = ELLI_DIRI;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "iter_err")) {
        par->i_err = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "max_iter")) {
        par->max_iter = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "elli_higen")) {
      par->method_itype = ELLI_HIGEN;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "distance")) {
         for (int i = 0; i < 4; i++) {
           par->distance[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "iter_err")) {
        par->i_err = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "max_iter")) {
        par->max_iter = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "parabolic")) {
      par->method_itype = PARABOLIC;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "o2i")) {
        par->o2i = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "hyperbolic")) {
      par->method_itype = HYPERBOLIC;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "epsilon")) {
        par->epsilon = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "bdry_type")) {
        par->bdry_itype = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "o2i")) {
        par->o2i = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "step_input_file")) {
        sprintf(par->step_input_file, "%s", thirditem->valuestring);
      }
    }
  }

  return ierr;
}


int
par_print(par_t *par)
{    
  int ierr = 0;

  fprintf(stdout,"input geometry file is \n %s\n",par->geometry_input_file);
  fprintf(stdout,"export grid dir is \n %s\n",par->grid_export_dir);
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
  if (par->check_ratio == 1) {
    fprintf(stdout, "------- check grid ratio-------\n");
  }
  if (par->check_step_x == 1) {
    fprintf(stdout, "------- check grid step x direction-------\n");
  }
  if (par->check_step_z == 1) {
    fprintf(stdout, "------- check grid step z direction-------\n");
  }
  if (par->check_smooth_x == 1) {
    fprintf(stdout, "------- check grid smooth x direction-------\n");
  }
  if (par->check_smooth_z == 1) {
    fprintf(stdout, "------- check grid smooth z direction-------\n");
  }

  if (par->flag_sample_x == 1) {
    fprintf(stdout,"sample grid x direction factor is %f\n",par->sample_factor_x);
  }
  if (par->flag_sample_z == 1) {
    fprintf(stdout,"sample grid z direction factor is %f\n",par->sample_factor_z);
  }
  if(par->flag_strech_x == 1) {
    fprintf(stdout, "------- strech x and strech coef is %f-------\n",par->strech_x_coef);
  }
  if(par->flag_strech_z == 1) {
    fprintf(stdout, "------- strech z and strech coef is %f-------\n",par->strech_z_coef);
  }

  fprintf(stdout, "------- grid generate method-------\n");
  if(par->method_itype == TFI) {
    fprintf(stdout, "grid generate method is linear TFI\n");
  }
  if(par->method_itype == HERMITE) {
    fprintf(stdout, "grid generate method is unidirection hermite\n");
    fprintf(stdout, "hermite coef is %f\n", par->coef);
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
  }
  if(par->method_itype == ELLI_DIRI) {
    fprintf(stdout, "grid generate method is elliptic_dirichlet\n");
    fprintf(stdout, "elli_diri coef is %f\n", par->coef);
    fprintf(stdout, "max_iteration is %d\n", par->max_iter);
    fprintf(stdout, "iter_error is %f\n", par->i_err);
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
  }
  if(par->method_itype == ELLI_HIGEN) {
    fprintf(stdout, "grid generate method is elliptic_hilgenstock\n");
    fprintf(stdout, "elli_higen coef is %f\n", par->coef);
    fprintf(stdout, "max_iteration is %d\n", par->max_iter);
    fprintf(stdout, "iter_error is %f\n", par->i_err);
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
    fprintf(stdout, "expect distance x1 bdry is %f\n",par->distance[0]);
    fprintf(stdout, "expect distance x2 bdry is %f\n",par->distance[1]);
    fprintf(stdout, "expect distance z1 bdry is %f\n",par->distance[2]);
    fprintf(stdout, "expect distance z2 bdry is %f\n",par->distance[3]);
  }
  if(par->method_itype == PARABOLIC) {
    fprintf(stdout, "grid generate method is parabolic\n");
    fprintf(stdout, "parabolic coef is %f\n", par->coef);
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
    if(par->o2i == 1)
    {
      fprintf(stdout, "outer(bdry_2) to inner(bdry_1)\n");
    } else {
      fprintf(stdout, "inner(bdry_1) to outer(bdry_2)\n");
    }
  }

  if(par->method_itype == HYPERBOLIC) {
    fprintf(stdout, "grid generate method is hyperbolic\n");
    fprintf(stdout, "hyperbolic coef is %f\n", par->coef);
    if(par->bdry_itype == 1) {
      fprintf(stdout, "boundary type is floating boundary\n");
    }
    if(par->bdry_itype == 2) {
      fprintf(stdout, "boundary type is cartesian boundary\n");
    }
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
    if(par->o2i == 1)
    {
      fprintf(stdout, "outer(bdry_2) to inner(bdry_1)\n");
    } else {
      fprintf(stdout, "inner(bdry_1) to outer(bdry_2)\n");
    }
    fprintf(stdout, "step file is  %s\n",par->step_input_file);
  }

  return ierr;
}
