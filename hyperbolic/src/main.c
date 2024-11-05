#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "par_t.h"
#include "constants.h"
#include "gd_t.h"
#include "io_funcs.h"
#include "quality_check.h"
#include "hyperbolic.h"

int main(int argc, char** argv)
{
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

  //-------------------------------------------------------------------------------
  // get commond-line argument
  //-------------------------------------------------------------------------------

  // argc checking
  if (argc < 2) {
    fprintf(stdout,"usage: main_grid_3d <par_file>\n");
    exit(1);
  }

  par_fname = argv[1];

  fprintf(stdout,"par file =  %s\n", par_fname); fflush(stdout);

  // read par
  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_read_from_file(par_fname, par);

  par_print(par);

  // generate grid 
  gd_t *gdcurv = (gd_t *) malloc(sizeof(gd_t));

  grid_init_set_hyper(gdcurv,par);

  time_t t_start = time(NULL);

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

  time_t t_end = time(NULL);

  fprintf(stdout,"\n************************************\n");
  fprintf(stdout,"grid generate running time is :%f s \n", difftime(t_end,t_start));
  fprintf(stdout,"************************************\n \n");

  fprintf(stdout,"export coord to file ... \n");
  gd_curv_coord_export(gdcurv,par);

  //  cal min dist 
  int indx_i, indx_j, indx_k;
  float dL_min;
  cal_min_dist(gdcurv, &indx_i, &indx_j, &indx_k, &dL_min);
  fprintf(stdout,"indx is (%d,%d,%d),dL_min_global is %f\n",
          indx_i, indx_j, indx_k, dL_min);
  fflush(stdout);

  // grid quality check and export quality data
  io_quality_t *io_quality = (io_quality_t *) malloc(sizeof(io_quality_t));
  if(par->grid_check == 1)
  {
    fprintf(stdout,"****************************************************** \n");
    fprintf(stdout,"***** grid quality check and export quality data ***** \n");
    fprintf(stdout,"****************************************************** \n");
    init_io_quality(io_quality,gdcurv);
    grid_quality_check(io_quality,gdcurv,par);
  }
  return 0;
}
