#ifndef SOLVER_H
#define SOLVER_H

/*************************************************
 * function prototype
 *************************************************/
int
thomas(int n, float *a, float *b, float *c, float *d_x, 
       float *d_y, float *d_z, float *u_x, float *u_y, float *u_z);

int
thomas_block(int n, double *a, double *b, double *c, double *d,
             double *x, double *D, double *y);

#endif
