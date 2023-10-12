#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "par_t.h"
#include "gd_t.h"
#include "io_funcs.h"
#include "algebra.h"
#include "constants.h"

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

  gd_t *gdcurv = (gd_t *) malloc(sizeof(gd_t));
  // for grid sample space
  gd_t *gdcurv_new = (gd_t *) malloc(sizeof(gd_t));

  read_import_coord(gdcurv,par);
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

  if(par->flag_sample == 1)
  {
    fprintf(stdout,"\n\n******** sample grid ******* \n");
    fprintf(stdout,"export coord to file ... \n");
    fprintf(stdout,"\n\n\n!!!!!!!!! NOTE !!!!!!!!!!!! \n");
    fprintf(stdout,"after sample grid, each direction abs(pml) layers must \n");
    fprintf(stdout,"multiply by sample factor\n");

    fprintf(stdout,"\n ****old abs layers ****\n");
    fprintf(stdout,"old number of pml layers x1 is %d\n",par->number_of_pml_x1);
    fprintf(stdout,"old number of pml layers x2 is %d\n",par->number_of_pml_x2);
    fprintf(stdout,"old number of pml layers y1 is %d\n",par->number_of_pml_y1);
    fprintf(stdout,"old number of pml layers y2 is %d\n",par->number_of_pml_y2);
    fprintf(stdout,"old number of pml layers z1 is %d\n",par->number_of_pml_z1);
    fprintf(stdout,"old number of pml layers z2 is %d\n",par->number_of_pml_z2);
    fprintf(stdout,"\n ****after sample new abs layers ****\n");
    fprintf(stdout,"new number of pml layers x1 is %d\n",par->number_of_pml_x1*par->sample_factor_xi);
    fprintf(stdout,"new number of pml layers x2 is %d\n",par->number_of_pml_x2*par->sample_factor_xi);
    fprintf(stdout,"new number of pml layers y1 is %d\n",par->number_of_pml_y1*par->sample_factor_et);
    fprintf(stdout,"new number of pml layers y2 is %d\n",par->number_of_pml_y2*par->sample_factor_et);
    fprintf(stdout,"new number of pml layers z1 is %d\n",par->number_of_pml_z1*par->sample_factor_zt);
    fprintf(stdout,"new number of pml layers z2 is %d\n",par->number_of_pml_z2*par->sample_factor_zt);
    fprintf(stdout,"currently sampling, please wait");
    fflush(stdout);

    grid_sample(gdcurv_new,gdcurv,par->sample_factor_xi,par->sample_factor_et,par->sample_factor_zt);
    gd_curv_coord_export(gdcurv_new,par);
  } else {
    fprintf(stdout,"******* not sample grid ******* \n");
    fprintf(stdout,"export coord to file ... \n");
    fflush(stdout);
    gd_curv_coord_export(gdcurv,par);
  }

  return 0;
}
