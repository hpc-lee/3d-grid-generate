#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H

#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
hyper_gene(gd_t *gdcurv, par_t *par);

int
cal_smooth_coef(gd_t *gdcurv, float coef, int k, float *coef_e_xi, float *coef_e_et);

int 
cal_matrix(gd_t *gdcurv, int k, float *a_xi, float *b_xi, float *c_xi, 
           float *a_et, float *b_et, float *c_et, float *d_et, float *volume);

int
modify_smooth(gd_t *gdcurv, int k, float *coef_e_xi, float *coef_e_et, 
              float *a_xi, float *b_xi, float *c_xi, float *a_et, 
              float *b_et, float *c_et, float *d_et);

int
modify_bdry(float *a_xi, float *b_xi, float *c_xi, float *a_et, 
            float *b_et, float *c_et, float *d_et, int nx, int ny,
            float epsilon_x, int bdry_x_itype,
            float epsilon_y, int bdry_y_itype);

int
assign_coords(gd_t *gdcurv, float *xyz, int k, float epsilon_x, 
              int bdry_x_itype, float epsilon_y, int bdry_y_itype);

int
solve_et_block(float *a_et, float *b_et, float *c_et,
               float *d_et, float *g_xi, int nx, int ny);

int
solve_xi_block(float *a_xi, float *b_xi, float *c_xi,
               float *g_xi, float *xyz, int nx, int ny);

 #endif
