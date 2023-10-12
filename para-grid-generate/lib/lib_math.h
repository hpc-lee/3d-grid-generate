#ifndef FDLIB_MATH_H
#define FDLIB_MATH_H

int
mat_invert3x3(double matrix[][3]);

int
mat_mul3x3(double A[][3], double B[][3], double C[][3]);

int
cross_product(double *A, double *B, double *C);

double
dot_product(double *A, double *B);

int 
mat_mul3x1(double A[][3], double *B, double *C);

int
mat_add3x3(double A[][3], double B[][3], double C[][3]);

int
vec_add3x1(double *A, double *B, double *C);

int
vec_sub3x1(double *A, double *B, double *C);

int
mat_sub3x3(double A[][3], double B[][3], double C[][3]);

int
mat_copy3x3(double A[][3], double B[][3]);

int
mat_iden3x3(double A[][3]);
#endif
