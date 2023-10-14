#ifndef SOLVER_H
#define SOLVER_H

/*************************************************
 * function prototype
 *************************************************/

int
update_SOR(float *x3d, float *y3d, float *z3d, float *x3d_tmp,
           float *y3d_tmp, float *z3d_tmp, int nx, int ny, int nz,
           float *P, float *Q, float *R, float w);

#endif
