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
#include "parabolic.h"

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

  // init MPI

  int myid, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &mpi_size);


  if (myid==0 && verbose>0) fprintf(stdout,"comm=%d, size=%d\n", comm, mpi_size); 
  if (myid==0 && verbose>0) fprintf(stdout,"par file =  %s\n", par_fname); 

  // read par
  par_t *par = (par_t *) malloc(sizeof(par_t));
  par_mpi_get(par_fname, myid, comm, par, verbose);
  if (myid==0 && verbose>0) par_print(par);

  // set mpi
  mympi_t *mympi = (mympi_t *) malloc(sizeof(mympi_t));
  if (myid==0 && verbose>0) fprintf(stdout,"set mpi topo ...\n"); 
  mympi_set(mympi,
            par->number_of_mpiprocs,
            comm,
            myid, verbose);

  gd_t *gdcurv = (gd_t *) malloc(sizeof(gd_t));
  // set gdinfo
  gd_info_set(gdcurv, mympi, par, verbose);

  gd_info_print(gdcurv, mympi);
  init_gdcurv(gdcurv);

  // set str in blk
  set_output_dir(gdcurv, mympi,
                 par->grid_export_dir,
                 verbose);
  
  // read bdry and init iter grid 
  bdry_t *bdry = (bdry_t *) malloc(sizeof(bdry_t));
  init_bdry(bdry,gdcurv);
  read_bdry(myid,bdry,par->geometry_input_file);
  assign_bdry_coord(gdcurv, bdry, mympi);

  grid_mesg_init(mympi, gdcurv);
  time_t t_start = time(NULL);
  para_gene(gdcurv, mympi, bdry, par);
  time_t t_end = time(NULL);
  if(myid == 0)
  {
    fprintf(stdout,"\n************************************\n");
    fprintf(stdout,"grid generate running time is :%f s \n", difftime(t_end,t_start));
    fprintf(stdout,"************************************\n \n");

    fprintf(stdout,"export coord to file ... \n");
  }
  gd_curv_coord_export(gdcurv,mympi);

  // grid quality check and export quality data
  io_quality_t *io_quality = (io_quality_t *) malloc(sizeof(io_quality_t));
  if(par->grid_check == 1)
  {
    if(myid == 0)
    {
      fprintf(stdout,"****************************************************** \n");
      fprintf(stdout,"***** grid quality check and export quality data ***** \n");
      fprintf(stdout,"****************************************************** \n");
    }
    init_io_quality(io_quality,gdcurv);
    grid_quality_check(io_quality,gdcurv,par,mympi);
  }

  MPI_Finalize();

  return 0;
}
