#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mympi_t.h"
#include "par_t.h"

//
// set grid size
//
int
mympi_set(mympi_t *mympi,
          int number_of_mpiprocs,
          MPI_Comm comm,
          int myid)
{
  int ierr = 0;
  
  mympi->nproc = number_of_mpiprocs;

  mympi->myid = myid;
  mympi->comm = comm;

  // mpi 3D topo 
  int pdims[1]   = {mympi->nproc};
  int periods[1] = {0};

  // create Cartesian topology
  MPI_Cart_create(comm, 1, pdims, periods, 0, &(mympi->topocomm));

  // get my local x,y coordinates
  MPI_Cart_coords(mympi->topocomm, myid, 1, mympi->topoid);

  // neighour
  MPI_Cart_shift(mympi->topocomm, 0, 1, &(mympi->neighid[0]), &(mympi->neighid[1]));

  return ierr;
}

