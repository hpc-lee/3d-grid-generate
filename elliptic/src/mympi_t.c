#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mympi_t.h"

//
// set grid size
//
int
mympi_set(mympi_t *mympi,
          int number_of_mpiprocs_x,
          int number_of_mpiprocs_y,
          int number_of_mpiprocs_z,
          MPI_Comm comm,
          int myid)
{
  int ierr = 0;

  mympi->nprocx = number_of_mpiprocs_x;
  mympi->nprocy = number_of_mpiprocs_y;
  mympi->nprocz = number_of_mpiprocs_z;

  mympi->myid = myid;
  mympi->comm = comm;

  // mpi topo, only consider 2d topo
  int pdims[3]   = {number_of_mpiprocs_x, number_of_mpiprocs_y, number_of_mpiprocs_z};
  int periods[3] = {0,0,0};

  // create Cartesian topology
  MPI_Cart_create(comm, 3, pdims, periods, 0, &(mympi->topocomm));

  // get my local x,y coordinates
  MPI_Cart_coords(mympi->topocomm, myid, 3, mympi->topoid);

  // neighour
  MPI_Cart_shift(mympi->topocomm, 0, 1, &(mympi->neighid[0]), &(mympi->neighid[1]));
  MPI_Cart_shift(mympi->topocomm, 1, 1, &(mympi->neighid[2]), &(mympi->neighid[3]));
  MPI_Cart_shift(mympi->topocomm, 2, 1, &(mympi->neighid[4]), &(mympi->neighid[5]));

  return ierr;
}

