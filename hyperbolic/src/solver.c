#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"
#include "lib_math.h"
#include "lib_mem.h"
#include "constants.h"

/*
  solve block tridiagonal linear system equaltion
  using thomas method
  Due to the sparse coefficient matrix of the tridiagonal equation,
  the computational complexity is proportional to n,
  rather than the n^3 of Gaussian elimination 

  [b1 c1        ]
  |a2 b2 c2     |
  |   a3 b3 c3  |
  [             ]

  a = square matrix element of lower diagonal
  b = square matrix element of main diagonal
  c = square matrix element of up diagonal
  d = right hand item, each element size is n*1

    [G1 0  0       ]
    |a2 G2 0       |
    |   a3 G3      |
  L=|              |
    |              |
    |              |
    [         an Gn]
                    
    [I1 D1 0            ]
    |0  I2 D2           |
    |   0  I3           |
  U=|                   |
    |                   |
    |          In-1 Dn-1|
    [          0    In  ]

  LU*x = d
  Ly=d
  Ux=y
*/
int
thomas_block(int n, double *a, double *b, double *c, double *d,
             double *x, double *D, double *y)
{
  double mat_a[3][3], mat_b[3][3], mat_c[3][3], vec_d[3];
  double mat_G[3][3], mat_D[3][3], vec_y[3], vec_x[3];
  double mat1[3][3],vec[3],vec1[3],vec2[3],vec3[3];
  double jac;
  size_t iptr1,iptr2,iptr3,iptr4;

  // i=0
  for(int ii=0; ii<3; ii++) {
    for(int jj=0; jj<3; jj++) {
      mat_G[ii][jj] = b[ii*3+jj];
      mat_c[ii][jj] = c[ii*3+jj];
    }
    vec_d[ii] = d[ii];
  }

  vec1[0]=mat_G[0][0]; vec1[1]=mat_G[0][1];vec1[2]=mat_G[0][2];
  vec2[0]=mat_G[1][0]; vec2[1]=mat_G[1][1];vec2[2]=mat_G[1][2];
  vec3[0]=mat_G[2][0]; vec3[1]=mat_G[2][1];vec3[2]=mat_G[2][2];
  cross_product(vec1,vec2,vec);
  jac = dot_product(vec,vec3);
  if(fabs(jac)<0.000001)
  {
    fprintf(stderr,"b{1} is singular, stop calculate");
    exit(1);
  } 

  mat_invert3x3(mat_G);
  mat_mul3x3(mat_G,mat_c,mat_D);
  mat_mul3x1(mat_G,vec_d,vec_y);
  for(int ii=0; ii<3; ii++) { 
    for(int jj=0; jj<3; jj++) {
      D[ii*3+jj] = mat_D[ii][jj];
    }
    y[ii] = vec_y[ii];
  }

  // i=1~n-1
  for(int i=1; i<n; i++)
  {
    iptr1 = i*3*3;
    iptr2 = i*3;
    iptr3 = (i-1)*3*3;
    iptr4 = (i-1)*3;
    for(int ii=0; ii<3; ii++) { 
      for(int jj=0; jj<3; jj++) {
        mat_a[ii][jj] = a[iptr1+ii*3+jj];
        mat_b[ii][jj] = b[iptr1+ii*3+jj];
        mat_c[ii][jj] = c[iptr1+ii*3+jj];
        mat_D[ii][jj] = D[iptr3+ii*3+jj];
      }
      vec_d[ii] = d[iptr2+ii];
      vec_y[ii] = y[iptr4+ii];
    }
    mat_mul3x3(mat_a,mat_D,mat1);  // a(i)*D(i-1)
    mat_sub3x3(mat_b,mat1,mat_G);  // G(i) = b(i)-a(i)*D(i-1)
    mat_invert3x3(mat_G);          // inv(G(i))
    mat_mul3x3(mat_G,mat_c,mat_D); // D(i) = inv(G(i))*c(i)
    mat_mul3x1(mat_a,vec_y,vec1);  // a(i) * y(i-1)
    vec_sub3x1(vec_d,vec1,vec2);   // d(i) - a(i) * y(i-1)
    mat_mul3x1(mat_G,vec2,vec_y);  // y(i) = inv(G(i))*(d(i) - a(i) * y(i-1))

    for(int ii=0; ii<3; ii++) { 
      for(int jj=0; jj<3; jj++) {
        D[iptr1+ii*3+jj] = mat_D[ii][jj]; 
      }
      y[iptr2+ii] = vec_y[ii];
    }
  }
  
  // i=n-1
  iptr2 = (n-1)*3;
  for(int ii=0; ii<3; ii++) { 
    x[iptr2+ii] = y[iptr2+ii];
  }

  for(int i=n-2; i>=0; i--)
  {
    iptr1 = i*3*3;
    iptr2 = i*3;
    iptr4 = (i+1)*3;
    
    for(int ii=0; ii<3; ii++) { 
      for(int jj=0; jj<3; jj++) {
        mat_D[ii][jj] = D[iptr1+ii*3+jj]; 
      }
      vec_x[ii] = x[iptr4+ii];
    }
    mat_mul3x1(mat_D,vec_x,vec1);  // D(i)*x(i+1)
    for(int ii=0; ii<3; ii++) 
    { 
      x[iptr2+ii] = y[iptr2+ii] - vec1[ii]; // u(i) = y(i) - D(i)*x(i+1) 
    }
  }
  
  return 0;
}

