#ifndef MY_MPI_H
#define MY_MPI_H

#include <mpi.h>

#include "constants.h"

/*******************************************************************************
 * structure
 ******************************************************************************/

typedef struct {

  int       nprocx;
  int       nprocy;
  int       nprocz;

  int       myid;
  MPI_Comm  comm;

  int    topoid[3];
  int    neighid[6];
  MPI_Comm    topocomm;

  MPI_Request *s_reqs;
  MPI_Request *r_reqs;

  int siz_sbuff_y1;
  int siz_sbuff_y2;

  int siz_rbuff_y1;
  int siz_rbuff_y2;

  int siz_sbuff;
  int siz_rbuff;

  float *sbuff;
  float *rbuff;

} mympi_t;

/*******************************************************************************
 * function prototype
 ******************************************************************************/

int
mympi_set(mympi_t *mympi,
          int number_of_mpiprocs,
          MPI_Comm comm, 
          int myid);
int
modify_neighid(mympi_t *mympi);

#endif
