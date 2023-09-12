#ifndef FDLIB_MATH_H
#define FDLIB_MATH_H

void
mat_invert2x2(double matrix[2][2]);

void
mat_mul2x2(double A[][2], double B[][2], double C[][2]);

void
mat_mul2x1(double A[][2], double B[2], double C[2]);

void
mat_add2x2(double A[][2], double B[][2], double C[][2]);

void
vec_add2x1(double A[2], double B[2], double C[2]);

void
mat_sub2x2(double A[][2], double B[][2], double C[][2]);

void
vec_sub2x1(double A[2], double B[2], double C[2]);

void
mat_copy2x2(double A[][2], double B[][2]);

void
mat_iden2x2(double A[][2]);

#endif
