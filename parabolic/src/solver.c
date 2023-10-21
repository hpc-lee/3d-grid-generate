#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"
#include "lib_math.h"
#include "lib_mem.h"
#include "constants.h"

// solve tridiagonal linear system equaltion
// using thomas method
// Due to the sparse coefficient matrix of the tridiagonal equation,
// the computational complexity is proportional to n,
// rather than the n^3 of Gaussian elimination
//
// [b1 c1        ]
// |a2 b2 c2     |
// |   a3 b3 c3  |
// [             ]
//
// a = element of lower diagonal
// b = element of main diagonal
// c = element of up diagonal
// d = right hand item
// u is unknow vector
int
thomas(int n, double *a, double *b, double *c, double *d_x, 
       double *d_y, double *d_z, double *u_x, double *u_y, double *u_z)
{
  double factor;
  for(int i=1; i<n; i++)
  {
    factor = a[i]/b[i-1];
    b[i] = b[i] - factor*c[i-1];
    d_x[i] = d_x[i] - factor*d_x[i-1];
    d_y[i] = d_y[i] - factor*d_y[i-1];
    d_z[i] = d_z[i] - factor*d_z[i-1];
  }

  u_x[n-1] = d_x[n-1]/b[n-1];
  u_y[n-1] = d_y[n-1]/b[n-1];
  u_z[n-1] = d_z[n-1]/b[n-1];
  for(int i=n-2; i>=0; i--)
  {
    u_x[i] = (d_x[i]-c[i]*u_x[i+1])/b[i];
    u_y[i] = (d_y[i]-c[i]*u_y[i+1])/b[i];
    u_z[i] = (d_z[i]-c[i]*u_z[i+1])/b[i];
  }

  return 0;
}

