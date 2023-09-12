#include <stdio.h>
#include <math.h>

#include "lib_math.h"

void
mat_invert2x2(double matrix[2][2])
{
  for (int k=0; k<2; k++)
  {
     double con = matrix[k][k];
     matrix[k][k] = 1.0;

     for (int i=0; i<2; i++) {
       matrix[k][i] = matrix[k][i]/con;
     }

     for (int i=0; i<2; i++)
     {
        if (i!=k) {
           con = matrix[i][k];
           matrix[i][k] = 0.0;
           for (int j=0; j<2; j++) {
             matrix[i][j] = matrix[i][j] - matrix[k][j] * con;
           }
        }
     }
  }

  return;
}

void
mat_mul2x2(double A[][2], double B[][2], double C[][2])
{
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++){
      C[i][j] = 0.0;
      for (int k=0; k<2; k++)
        C[i][j] += A[i][k] * B[k][j];
    }

  return;
}

void
mat_mul2x1(double A[][2], double B[2], double C[2])
{
  C[0] = A[0][0]*B[0] + A[0][1]*B[1];
  C[1] = A[1][0]*B[0] + A[1][1]*B[1];
  return;
}

void
mat_add2x2(double A[][2], double B[][2], double C[][2])
{
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
        C[i][j] = A[i][j] + B[i][j];
    }
  }
  return;
}

void
vec_add2x1(double A[2], double B[2], double C[2])
{
  for (int i=0; i<2; i++){
    C[i] = A[i] + B[i];
  }
  return;
}

void
vec_sub2x1(double A[2], double B[2], double C[2])
{
  for (int i=0; i<2; i++){
    C[i] = A[i] - B[i];
  }
  return;
}

void
mat_sub2x2(double A[][2], double B[][2], double C[][2])
{
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
        C[i][j] = A[i][j] - B[i][j];
    }
  }
  return;
}

void
mat_copy2x2(double A[][2], double B[][2])
{
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
        B[i][j] = A[i][j];
    }
  }
  return;
}

void
mat_iden2x2(double A[][2])
{
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      A[i][j] = 0.0;
      if(i==j) {
        A[i][j] = 1;
      }
    }
  }

  return;
}
