#ifndef SOLVER_H
#define SOLVER_H

/*************************************************
 * function prototype
 *************************************************/
int
thomas(int n, float *a, float *b, float *c, float *d_x, 
       float *d_y, float *d_z, float *u_x, float *u_y, float *u_z);

int
update_SOR(float *x3d, float *y3d, float *z3d, float *x3d_tmp,
           float *y3d_tmp, float *z3d_tmp, int nx, int ny, int nz,
           float *P, float *Q, float *R, float w);

int
thomas_block(int n, double *a, double *b, double *c, double *d,
             double *x, double *D, double *y);

#endif
