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
          int myid, int verbose)
{
  int ierr = 0;
  
  mympi->nprocx = 1;
  mympi->nprocy = number_of_mpiprocs;
  mympi->nprocz = 1;

  mympi->myid = myid;
  mympi->comm = comm;

  // mpi 3D topo 
  int pdims[3]   = {mympi->nprocx, mympi->nprocy, mympi->nprocz};
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

// when direction is y
// need swap neighid 
int
modify_neighid(mympi_t *mympi)
{
  int tmp_neighid_4 = mympi->neighid[4];
  int tmp_neighid_5 = mympi->neighid[5];
  mympi->neighid[4] = mympi->neighid[2];
  mympi->neighid[5] = mympi->neighid[3];
  mympi->neighid[2] = tmp_neighid_4;
  mympi->neighid[3] = tmp_neighid_5;

  int tmp_topodid = mympi->topoid[2];
  mympi->topoid[2] = mympi->topoid[1];
  mympi->topoid[1] = tmp_topodid;

  return 0;
}

