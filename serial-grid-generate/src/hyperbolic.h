#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H

#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
hyper_gene(gd_t *gdcurv, par_t *par);

int
cal_smooth_coef(float coef, float *x3d, float *y3d, float *z3d, int nx,
                int ny, int nz, int k, float *coef_e_xi, float *coef_e_et);

int 
cal_matrix(float *x3d,float *y3d,float *z3d, int nx, int ny, int k,
           float *step, double *a_xi, double *b_xi, double *c_xi, double *a_et,
           double *b_et, double *c_et, double *d_et, float *volume);

int
modify_smooth(float *x3d, float *y3d, float *z3d, int nx, int ny, int k,
              float *coef_e_xi, float *coef_e_et, double *a_xi, double *b_xi,
              double *c_xi, double *a_et, double *b_et, 
              double *c_et, double *d_et);

int
modify_bdry(double *a_xi, double *b_xi, double *c_xi, double *a_et, 
            double *b_et, double *c_et, double *d_et, int nx, int ny,
            float epsilon_x, int bdry_x_itype,
            float epsilon_y, int bdry_y_itype);

int
assign_coords(double *xyz, float *x3d, float *y3d, float *z3d,
              int nx, int ny,int k, float epsilon_x, int bdry_x_itype,
              float epsilon_y, int bdry_y_itype);

int
solve_et_block(double *a_et, double *b_et, double *c_et,
               double *d_et, double *g_xi, int nx, int ny);

int
solve_xi_block(double *a_xi, double *b_xi, double *c_xi,
               double *g_xi, double *xyz, int nx, int ny);

 #endif
