#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "par_t.h"
#include "constants.h"
#include "gd_t.h"
#include "algebra.h"
#include "io_funcs.h"
#include "quality_check.h"
#include "parabolic.h"
#include "elliptic.h"
#include "hyperbolic.h"

int main(int argc, char** argv)
{
  int verbose;
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

  //-------------------------------------------------------------------------------
  // get commond-line argument
  //-------------------------------------------------------------------------------

  // argc checking
  if (argc < 3) {
    fprintf(stdout,"usage: main_grid_3d <par_file> <opt: verbose>\n");
    exit(1);
  }

  par_fname = argv[1];

  if (argc >= 3) {
    verbose = atoi(argv[2]); // verbose number
    fprintf(stdout,"verbose=%d\n", verbose); fflush(stdout);
  }

  fprintf(stdout,"par file =  %s\n", par_fname); fflush(stdout);

  // read par
  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_read_from_file(par_fname, par, verbose);

  if (verbose>0) par_print(par);

  // generate grid 
  gd_t *gdcurv = (gd_t *) malloc(sizeof(gd_t));
  // for grid sample space
  gd_t *gdcurv_new = (gd_t *) malloc(sizeof(gd_t));
  switch(par->method_itype)
  {
    case TFI : {

      grid_init_set(gdcurv,par->geometry_input_file);
      // grid method
      linear_tfi(gdcurv);

      break;
    }

    case ELLI_DIRI : {

      grid_init_set(gdcurv,par->geometry_input_file);
      // linear tfi generate init iter grid
      linear_tfi(gdcurv);
      diri_gene(gdcurv,par);

      break;
    }
    case ELLI_HIGEN : {

      grid_init_set(gdcurv,par->geometry_input_file);
      // linear tfi generate init iter grid
      linear_tfi(gdcurv);
      higen_gene(gdcurv,par);

      break;
    }
    case PARABOLIC : {

      grid_init_set(gdcurv,par->geometry_input_file);
      // before grid generate
      if(par->dire_itype == X_DIRE)
      {
        permute_coord_x(gdcurv);
      }
      if(par->dire_itype == Y_DIRE)
      {
        permute_coord_y(gdcurv);
      }

      para_gene(gdcurv,par);

      // after grid generate
      if(par->dire_itype == X_DIRE)
      {
        permute_coord_x(gdcurv);
      }
      if(par->dire_itype == Y_DIRE)
      {
        permute_coord_y(gdcurv);
      }

      break;
    }
    case HYPERBOLIC : {

      grid_init_set_hyper(gdcurv,par);
      hyper_gene(gdcurv,par);

      // after grid generate
      if(par->dire_itype == X_DIRE)
      {
        permute_coord_x(gdcurv);
      }
      if(par->dire_itype == Y_DIRE)
      {
        permute_coord_y(gdcurv);
      }

      break;
    }
  }
  // strech xi-direction grid
  if(par->flag_strech_xi == 1)
  {
    xi_arc_strech(gdcurv,par->strech_xi_coef);
  }

  // strech et-direction grid
  if(par->flag_strech_et == 1)
  {
    et_arc_strech(gdcurv,par->strech_et_coef);
  }

  // strech zt-direction grid
  if(par->flag_strech_zt == 1)
  {
    zt_arc_strech(gdcurv,par->strech_zt_coef);
  }

  if(par->flag_sample_xi == 1 || par->flag_sample_et == 1 || par->flag_sample_zt == 1)
  {
    grid_sample(gdcurv_new,gdcurv,par->sample_factor_xi,par->sample_factor_et,par->sample_factor_zt);
    fprintf(stdout,"******** sample grid ******* \n");
    fprintf(stdout,"export coord to file ... \n");
    gd_curv_coord_export(gdcurv_new,par->grid_export_dir);
  } else {
    fprintf(stdout,"******* not sample grid ******* \n");
    fprintf(stdout,"export coord to file ... \n");
    gd_curv_coord_export(gdcurv,par->grid_export_dir);
  }

  // grid quality check and export quality data
  io_quality_t *io_quality = (io_quality_t *) malloc(sizeof(io_quality_t));
  if(par->grid_check == 1)
  {
    fprintf(stdout,"****************************************************** \n");
    fprintf(stdout,"***** grid quality check and export quality data ***** \n");
    fprintf(stdout,"****************************************************** \n");
    if(par->flag_sample_xi == 1 || par->flag_sample_et == 1 || par->flag_sample_zt == 1)
    {
      init_io_quality(io_quality,gdcurv_new);
      grid_quality_check(io_quality,gdcurv_new,par);
    } else {
      init_io_quality(io_quality,gdcurv);
      grid_quality_check(io_quality,gdcurv,par);
    }
  }
  return 0;
}