#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "gd_t.h"

typedef struct {

  // x1 bdry1
  double *dif_b1_x_k_loc; 
  double *dif_b1_y_k_loc; 
  double *dif_b1_z_k_loc; 
  double *dif_b1_x_k1_loc;
  double *dif_b1_y_k1_loc;
  double *dif_b1_z_k1_loc;

  double *dif_b1_x_k; 
  double *dif_b1_y_k; 
  double *dif_b1_z_k; 
  double *dif_b1_x_k1;
  double *dif_b1_y_k1;
  double *dif_b1_z_k1;

  // x2 bdry2
  double *dif_b2_x_k_loc; 
  double *dif_b2_y_k_loc; 
  double *dif_b2_z_k_loc; 
  double *dif_b2_x_k1_loc;
  double *dif_b2_y_k1_loc;
  double *dif_b2_z_k1_loc;

  double *dif_b2_x_k; 
  double *dif_b2_y_k; 
  double *dif_b2_z_k; 
  double *dif_b2_x_k1;
  double *dif_b2_y_k1;
  double *dif_b2_z_k1;

  // y1 bdry3
  double *dif_b3_x_k_loc; 
  double *dif_b3_y_k_loc; 
  double *dif_b3_z_k_loc; 
  double *dif_b3_x_k1_loc;
  double *dif_b3_y_k1_loc;
  double *dif_b3_z_k1_loc;

  double *dif_b3_x_k; 
  double *dif_b3_y_k; 
  double *dif_b3_z_k; 
  double *dif_b3_x_k1;
  double *dif_b3_y_k1;
  double *dif_b3_z_k1;

  // y2 bdry4
  double *dif_b4_x_k_loc; 
  double *dif_b4_y_k_loc; 
  double *dif_b4_z_k_loc; 
  double *dif_b4_x_k1_loc;
  double *dif_b4_y_k1_loc;
  double *dif_b4_z_k1_loc;

  double *dif_b4_x_k; 
  double *dif_b4_y_k; 
  double *dif_b4_z_k; 
  double *dif_b4_x_k1;
  double *dif_b4_y_k1;
  double *dif_b4_z_k1;

} bdry_effct_t;
/*************************************************
 * function prototype
 *************************************************/

int 
para_gene(gd_t *gdcurv, mympi_t *mympi,
          bdry_t *bdry, par_t *par);

int 
predict_point(gd_t *gdcurv, bdry_effct_t *bdry_effct, mympi_t *mympi,
              int k, int o2i, double coef, double *x1_len,
              double *x2_len, double *y1_len, double *y2_len);

int
update_point(gd_t *gdcurv, double *var_th, int k, double *coord); 

int
bdry_modify(gd_t *gdcurv, bdry_effct_t *bdry_effct, mympi_t *mympi, int k);

int
flip_bdry_z(double *x1, double *x2, double *y1, double *y2,
            int total_nx, int total_ny, int total_nz);

int
cal_bdry_arc_length(double *x1, double *x2, double *y1, double *y2, 
                    int total_nx, int total_ny, int total_nz, 
                    double *x1_len, double *x2_len, double *y1_len, double *y2_len);

int
exchange_coord(gd_t *gdcurv, mympi_t *mympi, int k, int num_of_s_reqs, int num_of_r_reqs);

int 
init_bdry_effct(bdry_effct_t *bdry_effct, gd_t *gdcurv);

#endif
