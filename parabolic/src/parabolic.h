#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "gd_t.h"

typedef struct {

  // x1 bdry1
  float *dif_b1_x_k_loc; 
  float *dif_b1_y_k_loc; 
  float *dif_b1_z_k_loc; 
  float *dif_b1_x_k1_loc;
  float *dif_b1_y_k1_loc;
  float *dif_b1_z_k1_loc;

  float *dif_b1_x_k; 
  float *dif_b1_y_k; 
  float *dif_b1_z_k; 
  float *dif_b1_x_k1;
  float *dif_b1_y_k1;
  float *dif_b1_z_k1;

  // x2 bdry2
  float *dif_b2_x_k_loc; 
  float *dif_b2_y_k_loc; 
  float *dif_b2_z_k_loc; 
  float *dif_b2_x_k1_loc;
  float *dif_b2_y_k1_loc;
  float *dif_b2_z_k1_loc;

  float *dif_b2_x_k; 
  float *dif_b2_y_k; 
  float *dif_b2_z_k; 
  float *dif_b2_x_k1;
  float *dif_b2_y_k1;
  float *dif_b2_z_k1;

  // y1 bdry3
  float *dif_b3_x_k_loc; 
  float *dif_b3_y_k_loc; 
  float *dif_b3_z_k_loc; 
  float *dif_b3_x_k1_loc;
  float *dif_b3_y_k1_loc;
  float *dif_b3_z_k1_loc;

  float *dif_b3_x_k; 
  float *dif_b3_y_k; 
  float *dif_b3_z_k; 
  float *dif_b3_x_k1;
  float *dif_b3_y_k1;
  float *dif_b3_z_k1;

  // y2 bdry4
  float *dif_b4_x_k_loc; 
  float *dif_b4_y_k_loc; 
  float *dif_b4_z_k_loc; 
  float *dif_b4_x_k1_loc;
  float *dif_b4_y_k1_loc;
  float *dif_b4_z_k1_loc;

  float *dif_b4_x_k; 
  float *dif_b4_y_k; 
  float *dif_b4_z_k; 
  float *dif_b4_x_k1;
  float *dif_b4_y_k1;
  float *dif_b4_z_k1;

} bdry_effct_t;
/*************************************************
 * function prototype
 *************************************************/

int 
para_gene(gd_t *gdcurv, mympi_t *mympi,
          bdry_t *bdry, par_t *par);

int 
predict_point(gd_t *gdcurv, bdry_effct_t *bdry_effct, mympi_t *mympi,
              int k, int o2i, float coef, float *x1_len,
              float *x2_len, float *y1_len, float *y2_len);

int
update_point(gd_t *gdcurv, float *var_th, int k, float *coord); 

int
bdry_modify(gd_t *gdcurv, bdry_effct_t *bdry_effct, mympi_t *mympi, int k);

int
flip_bdry_z(float *x1, float *x2, float *y1, float *y2,
            int total_nx, int total_ny, int total_nz);

int
cal_bdry_arc_length(float *x1, float *x2, float *y1, float *y2, 
                    int total_nx, int total_ny, int total_nz, 
                    float *x1_len, float *x2_len, float *y1_len, float *y2_len);

int
exchange_coord(gd_t *gdcurv, mympi_t *mympi, int k, int num_of_s_reqs, int num_of_r_reqs);

int 
init_bdry_effct(bdry_effct_t *bdry_effct, gd_t *gdcurv);

#endif
