#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "elliptic.h"
#include "constants.h"
#include "solver.h"
#include "gd_t.h"
#include "lib_mem.h"

int
init_src(src_t *src, gd_t *gdcurv)
{
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;
  
  // not include 2 bdry points
  size_t siz_bx = (total_ny-2)*(total_nz-2);
  size_t siz_by = (total_nx-2)*(total_nz-2);
  size_t siz_bz = (total_nx-2)*(total_ny-2);

  src->P_x1_loc = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x1");
  src->Q_x1_loc = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x1");
  src->R_x1_loc = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x1");
  src->P_x2_loc = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x2");
  src->Q_x2_loc = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x2");
  src->R_x2_loc = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x2");

  src->P_y1_loc = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y1");
  src->Q_y1_loc = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y1");
  src->R_y1_loc = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y1");
  src->P_y2_loc = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y2");
  src->Q_y2_loc = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y2");
  src->R_y2_loc = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y2");

  src->P_z1_loc = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z1");
  src->Q_z1_loc = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z1");
  src->R_z1_loc = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z1");
  src->P_z2_loc = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z2");
  src->Q_z2_loc = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z2");
  src->R_z2_loc = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z2");

  src->P_x1 = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x1");
  src->Q_x1 = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x1");
  src->R_x1 = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x1");
  src->P_x2 = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x2");
  src->Q_x2 = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x2");
  src->R_x2 = (float *)mem_calloc_1d_float(siz_bx, 0.0, "P_x2");

  src->P_y1 = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y1");
  src->Q_y1 = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y1");
  src->R_y1 = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y1");
  src->P_y2 = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y2");
  src->Q_y2 = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y2");
  src->R_y2 = (float *)mem_calloc_1d_float(siz_by, 0.0, "P_y2");

  src->P_z1 = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z1");
  src->Q_z1 = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z1");
  src->R_z1 = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z1");
  src->P_z2 = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z2");
  src->Q_z2 = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z2");
  src->R_z2 = (float *)mem_calloc_1d_float(siz_bz, 0.0, "P_z2");

  // include 2 ghost points
  src->P = (float *)mem_calloc_1d_float(nx*ny*nz, 0.0, "source P");
  src->Q = (float *)mem_calloc_1d_float(nx*ny*nz, 0.0, "source Q");
  src->R = (float *)mem_calloc_1d_float(nx*ny*nz, 0.0, "source R");
}

int
diri_gene(gd_t *gdcurv, par_t *par, mympi_t *mympi)
{
  float err = par->iter_err;
  int max_iter = par->max_iter;
  float *coef = par->coef;

  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t iptr, iptr1, iptr2, iptr3;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  float *x3d_tmp = gdcurv->x3d_tmp;
  float *y3d_tmp = gdcurv->y3d_tmp;
  float *z3d_tmp = gdcurv->z3d_tmp;

  MPI_Comm comm = mympi->comm;
  int myid = mympi->myid;
  int *neighid = mympi->neighid;

  int num_of_s_reqs = 6;
  int num_of_r_reqs = 6;
  
  src_t *src = (src_t *) malloc(sizeof(src_t));
  init_src(src,gdcurv);

  // before update grid. use init grid to
  // calculate ghost points and bdry g11,g22,g33
  float *p_x; // point_x_bdry
  float *p_y; // point_y_bdry
  float *p_z; // point_z_bdry
  float *g11_x;  // g11_x_bdry
  float *g22_y;  // g22_y_bdry
  float *g33_z;  // g33_z_bdry
  p_x = (float *)mem_calloc_1d_float(nz*ny*3*2, 0.0, "p_x");
  p_y = (float *)mem_calloc_1d_float(nz*nx*3*2, 0.0, "p_y");
  p_z = (float *)mem_calloc_1d_float(ny*nx*3*2, 0.0, "p_z");
  g11_x = (float *)mem_calloc_1d_float(nz*ny*2, 0.0, "g11_x");
  g22_y = (float *)mem_calloc_1d_float(nz*nx*2, 0.0, "g22_y");
  g33_z = (float *)mem_calloc_1d_float(ny*nx*2, 0.0, "g33_z");
  ghost_cal(x3d,y3d,z3d,nx,ny,nz,p_x,p_y,p_z,g11_x,g22_y,g33_z,neighid);

  set_src_diri(x3d,y3d,z3d,gdcurv,src,p_x,p_y,p_z,g11_x,g22_y,g33_z,mympi);
  interp_inner_source(src, gdcurv, coef);

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + j*siz_iy + i;
        x3d_tmp[iptr] = x3d[iptr];
        y3d_tmp[iptr] = y3d[iptr];
        z3d_tmp[iptr] = z3d[iptr];
      }
    }
  }

  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  float *ptr_x, *ptr_y, *ptr_z;
  
  int flag_true = 1;
  float resi, resj, resk;
  float max_resi, max_resj, max_resk;
  float dif, dif1, dif2, dif3, dif_x, dif_y, dif_z;
  while(flag_true)
  {
    // update solver
    update_SOR(x3d,y3d,z3d,x3d_tmp,y3d_tmp,z3d_tmp,nx,ny,nz,src->P,src->Q,src->R,w);
    Niter += 1;

    MPI_Startall(num_of_r_reqs, mympi->r_reqs);

    grid_pack_mesg(mympi,gdcurv,x3d_tmp,y3d_tmp,z3d_tmp);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs);

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs,MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs,MPI_STATUS_IGNORE);

    grid_unpack_mesg(mympi,gdcurv,x3d_tmp,y3d_tmp,z3d_tmp);

    max_resi = 0.0;
    max_resj = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + i+1;
          iptr2 = k*siz_iz + (j+1)*siz_iy  + i;
          iptr3 = (k+1)*siz_iz + j*siz_iy  + i;
          dif_x = x3d_tmp[iptr]-x3d[iptr];
          dif_y = y3d_tmp[iptr]-y3d[iptr];
          dif_z = z3d_tmp[iptr]-z3d[iptr];
          dif = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          dif_x = x3d[iptr1]-x3d[iptr];
          dif_y = y3d[iptr1]-y3d[iptr];
          dif_z = z3d[iptr1]-z3d[iptr];
          dif1 = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          dif_x = x3d[iptr2]-x3d[iptr];
          dif_y = y3d[iptr2]-y3d[iptr];
          dif_z = z3d[iptr2]-z3d[iptr];
          dif2 = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          dif_x = x3d[iptr3]-x3d[iptr];
          dif_y = y3d[iptr3]-y3d[iptr];
          dif_z = z3d[iptr3]-z3d[iptr];
          dif3 = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          resi = dif/dif1;
          resj = dif/dif2;
          resk = dif/dif3;
          max_resi = fmax(max_resi,resi);
          max_resj = fmax(max_resj,resj);
          max_resk = fmax(max_resk,resk);
        }
      }
    }

    float sendbuf = max_resi;
    MPI_Allreduce(&sendbuf, &max_resi, 1, MPI_REAL, MPI_MAX, comm);
    sendbuf = max_resj;
    MPI_Allreduce(&sendbuf, &max_resj, 1, MPI_REAL, MPI_MAX, comm);
    sendbuf = max_resk;
    MPI_Allreduce(&sendbuf, &max_resk, 1, MPI_REAL, MPI_MAX, comm);

    if(myid == 0)
    {
      fprintf(stdout,"number of iter is %d\n", Niter);
      fprintf(stdout,"max_resi is %f, max_resj is %f, max_resk is %f\n", max_resi, max_resj, max_resk);
    }

    // swap pointer, avoid coping
    ptr_x = x3d;
    ptr_y = y3d;
    ptr_z = z3d;
    x3d = x3d_tmp;
    y3d = y3d_tmp;
    z3d = z3d_tmp;
    x3d_tmp = ptr_x;  
    y3d_tmp = ptr_y;
    z3d_tmp = ptr_z;
    
    set_src_diri(x3d,y3d,z3d,gdcurv,src,p_x,p_y,p_z,g11_x,g22_y,g33_z,mympi);
    interp_inner_source(src, gdcurv, coef);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resj < err && max_resk < err)
    {
      flag_true = 0;
    }

  }
  
  gdcurv->v4d = x3d;
  gdcurv->x3d = x3d;
  gdcurv->y3d = y3d;
  gdcurv->z3d = z3d;

  gdcurv->v4d_tmp = x3d_tmp;

  free(gdcurv->v4d_tmp);  // free temp grid space
  free(p_x);
  free(p_y);
  free(p_z);
  free(g11_x);
  free(g22_y);
  free(g33_z);

  return 0;
}

int
set_src_diri(float *x3d, float *y3d, float *z3d, gd_t *gdcurv, 
             src_t *src, float *p_x, float *p_y, float *p_z,
             float *g11_x, float *g22_y, float *g33_z, 
             mympi_t *mympi)
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  int gni1 = gdcurv->gni1;
  int gnj1 = gdcurv->gnj1;
  int gnk1 = gdcurv->gnk1;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int gni,gnj,gnk;
  size_t size;

  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float temp_x,temp_y,temp_z; 
  float temp1,temp2,temp3; 
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xixi,y_xixi,z_xixi;
  float x_etet,y_etet,z_etet;
  float x_ztzt,y_ztzt,z_ztzt;
  float x_xizt,y_xizt,z_xizt;
  float x_xiet,y_xiet,z_xiet;
  float x_etzt,y_etzt,z_etzt;
  float g11,g22,g33,g12,g13,g23;
  float c1,c2,c3,c4;

  float *p_x1,*p_x2;
  float *p_y1,*p_y2;
  float *p_z1,*p_z2;
  float *g11_x1,*g11_x2;
  float *g22_y1,*g22_y2;
  float *g33_z1,*g33_z2;
  p_x1 = p_x;
  p_x2 = p_x + ny*nz*3;
  p_y1 = p_y;
  p_y2 = p_y + nx*nz*3;
  p_z1 = p_z;
  p_z2 = p_z + nx*ny*3;
  g11_x1 = g11_x;
  g11_x2 = g11_x + ny*nz;
  g22_y1 = g22_y;
  g22_y2 = g22_y + nx*nz;
  g33_z1 = g33_z;
  g33_z2 = g33_z + nx*ny;

  MPI_Comm topocomm = mympi->topocomm;
  int *neighid = mympi->neighid;

  float *P_x1_loc = src->P_x1_loc;
  float *Q_x1_loc = src->Q_x1_loc;
  float *R_x1_loc = src->R_x1_loc;
  float *P_x2_loc = src->P_x2_loc;
  float *Q_x2_loc = src->Q_x2_loc;
  float *R_x2_loc = src->R_x2_loc;

  float *P_y1_loc = src->P_y1_loc;
  float *Q_y1_loc = src->Q_y1_loc;
  float *R_y1_loc = src->R_y1_loc;
  float *P_y2_loc = src->P_y2_loc;
  float *Q_y2_loc = src->Q_y2_loc;
  float *R_y2_loc = src->R_y2_loc;
  
  float *P_z1_loc = src->P_z1_loc;
  float *Q_z1_loc = src->Q_z1_loc;
  float *R_z1_loc = src->R_z1_loc;
  float *P_z2_loc = src->P_z2_loc;
  float *Q_z2_loc = src->Q_z2_loc;
  float *R_z2_loc = src->R_z2_loc;

  float *P_x1 = src->P_x1;
  float *Q_x1 = src->Q_x1;
  float *R_x1 = src->R_x1;
  float *P_x2 = src->P_x2;
  float *Q_x2 = src->Q_x2;
  float *R_x2 = src->R_x2;

  float *P_y1 = src->P_y1;
  float *Q_y1 = src->Q_y1;
  float *R_y1 = src->R_y1;
  float *P_y2 = src->P_y2;
  float *Q_y2 = src->Q_y2;
  float *R_y2 = src->R_y2;
  
  float *P_z1 = src->P_z1;
  float *Q_z1 = src->Q_z1;
  float *R_z1 = src->R_z1;
  float *P_z2 = src->P_z2;
  float *Q_z2 = src->Q_z2;
  float *R_z2 = src->R_z2;
  // bdry x1 xi=0
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr  = k*siz_iz + j*siz_iy + 0;     //(0,j,k)
        iptr1 = k*siz_iz + (j+1)*siz_iy + 0; //(0,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + 0; //(0,j-1,k)
        iptr3 = (k+1)*siz_iz + j*siz_iy + 0; //(0,j,k+1)
        iptr4 = (k-1)*siz_iz + j*siz_iy + 0; //(0,j,k-1)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]); 
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]); 
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]); 
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]); 
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]); 
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]); 
        x_etet = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
        y_etet = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
        z_etet = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
        x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
        y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
        z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

        iptr1 = (k-1)*siz_iz + (j-1)*siz_iy + 0; //(0,j-1,k-1)
        iptr2 = (k-1)*siz_iz + (j+1)*siz_iy + 0; //(0,j+1,k-1)
        iptr3 = (k+1)*siz_iz + (j+1)*siz_iy + 0; //(0,j+1,k+1)
        iptr4 = (k+1)*siz_iz + (j-1)*siz_iy + 0; //(0,j-1,k+1)

        x_etzt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
        y_etzt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
        z_etzt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
        
        iptr1 = k*siz_iz + j*siz_iy + 1;     //(1,j,k)
        iptr2 = k*ny + j;
        size = ny*nz;
        x_xi = x3d[iptr1] - x3d[iptr];
        y_xi = y3d[iptr1] - y3d[iptr];
        z_xi = z3d[iptr1] - z3d[iptr];
        x_xixi = p_x1[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
        y_xixi = p_x1[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
        z_xixi = p_x1[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

        g22 = x_et*x_et + y_et*y_et + z_et*z_et;
        g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
        g23 = x_et*x_zt + y_et*y_zt + z_et*z_zt;

        c1 = 1/(g22*g33 - g23*g23); // 1/alpha1
        c2 = g23/g33;
        c3 = g23/g22;
        c4 = 1/g11_x1[iptr2];

        temp1 = -c4*(x_xi*x_xixi + y_xi*y_xixi + z_xi*z_xixi);
        temp2 = -c4*((x_et-c2*x_zt)*x_xixi + (y_et-c2*y_zt)*y_xixi + (z_et-c2*z_zt)*z_xixi);
        temp3 = -c4*((x_zt-c3*x_et)*x_xixi + (y_zt-c3*y_et)*y_xixi + (z_zt-c3*z_et)*z_xixi);
        temp_x = g33*x_etet + g22*x_ztzt - 2*g23*x_etzt;
        temp_y = g33*y_etet + g22*y_ztzt - 2*g23*y_etzt;
        temp_z = g33*z_etet + g22*z_ztzt - 2*g23*z_etzt;
        
        gnj = gnj1 + (j-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_ny-2) + gnj;
        P_x1_loc[iptr] = temp1 - c1*(x_xi*temp_x + y_xi*temp_y + z_xi*temp_z);

        Q_x1_loc[iptr] = temp2 - c1*((x_et-c2*x_zt)*temp_x + (y_et-c2*y_zt)*temp_y 
                       + (z_et-c2*z_zt)*temp_z);

        R_x1_loc[iptr] = temp3 - c1*((x_zt-c3*x_et)*temp_x + (y_zt-c3*y_et)*temp_y 
                       + (z_zt-c3*z_et)*temp_z);
      }
    }
  }

  // bdry x2 
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr = k*siz_iz + j*siz_iy + nx-1;
        iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1; //(nx-1,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1; //(nx-1,j-1,k)
        iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k+1)
        iptr4 = (k-1)*siz_iz + j*siz_iy + nx-1; //(nx-1,j,k-1)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
        x_etet = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
        y_etet = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
        z_etet = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
        x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
        y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
        z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

        iptr1 = (k-1)*siz_iz + (j-1)*siz_iy + nx-1; //(nx-1,j-1,k-1)
        iptr2 = (k-1)*siz_iz + (j+1)*siz_iy + nx-1; //(nx-1,j+1,k-1)
        iptr3 = (k+1)*siz_iz + (j+1)*siz_iy + nx-1; //(nx-1,j+1,k+1)
        iptr4 = (k+1)*siz_iz + (j-1)*siz_iy + nx-1; //(nx-1,j-1,k+1)

        x_etzt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
        y_etzt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
        z_etzt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
        
        iptr1 = k*siz_iz + j*siz_iy + nx-2;     //(nx-2,j,k)
        iptr2 = k*ny + j;
        size = ny*nz;
        x_xi = x3d[iptr] - x3d[iptr1];
        y_xi = y3d[iptr] - y3d[iptr1];
        z_xi = z3d[iptr] - z3d[iptr1];
        x_xixi = p_x2[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
        y_xixi = p_x2[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
        z_xixi = p_x2[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

        g22 = x_et*x_et + y_et*y_et + z_et*z_et;
        g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
        g23 = x_et*x_zt + y_et*y_zt + z_et*z_zt;

        c1 = 1/(g22*g33 - g23*g23); // 1/alpha1
        c2 = g23/g33;
        c3 = g23/g22;
        c4 = 1/g11_x2[iptr2];

        temp1 = -c4*(x_xi*x_xixi + y_xi*y_xixi + z_xi*z_xixi);
        temp2 = -c4*((x_et-c2*x_zt)*x_xixi + (y_et-c2*y_zt)*y_xixi + (z_et-c2*z_zt)*z_xixi);
        temp3 = -c4*((x_zt-c3*x_et)*x_xixi + (y_zt-c3*y_et)*y_xixi + (z_zt-c3*z_et)*z_xixi);
        temp_x = g33*x_etet + g22*x_ztzt - 2*g23*x_etzt;
        temp_y = g33*y_etet + g22*y_ztzt - 2*g23*y_etzt;
        temp_z = g33*z_etet + g22*z_ztzt - 2*g23*z_etzt;

        gnj = gnj1 + (j-1);
        gnk = gnk1 + (k-1);
        iptr  = gnk*(total_ny-2) + gnj;  

        P_x2_loc[iptr] = temp1 - c1*(x_xi*temp_x + y_xi*temp_y + z_xi*temp_z);

        Q_x2_loc[iptr] = temp2 - c1*((x_et-c2*x_zt)*temp_x + (y_et-c2*y_zt)*temp_y 
                       + (z_et-c2*z_zt)*temp_z);

        R_x2_loc[iptr] = temp3 - c1*((x_zt-c3*x_et)*temp_x + (y_zt-c3*y_et)*temp_y 
                       + (z_zt-c3*z_et)*temp_z);

      }
    }
  }

  // bdry y1 
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*siz_iz + 0*siz_iy + i;   //(i,0,k)
        iptr1 = k*siz_iz + 0*siz_iy + i+1;   //(i+1,0,k)
        iptr2 = k*siz_iz + 0*siz_iy + i-1;   //(i-1,0,k)
        iptr3 = (k+1)*siz_iz + 0*siz_iy + i; //(i,0,k+1)
        iptr4 = (k-1)*siz_iz + 0*siz_iy + i; //(i,0,k-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]); 
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]); 
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]); 
        x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
        y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
        z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
        x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
        y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
        z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

        iptr1 = (k-1)*siz_iz + 0*siz_iy + i-1; //(i-1,0,k-1)
        iptr2 = (k-1)*siz_iz + 0*siz_iy + i+1; //(i+1,0,k-1)
        iptr3 = (k+1)*siz_iz + 0*siz_iy + i+1; //(i+1,0,k+1)
        iptr4 = (k+1)*siz_iz + 0*siz_iy + i-1; //(i-1,0,k+1)

        x_xizt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
        y_xizt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
        z_xizt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
        
        iptr1 = k*siz_iz + 1*siz_iy + i;     //(i,1,k)
        iptr2 = k*nx + i;
        size = nx*nz;
        x_et = x3d[iptr1] - x3d[iptr];
        y_et = y3d[iptr1] - y3d[iptr];
        z_et = z3d[iptr1] - z3d[iptr];
        x_etet = p_y1[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
        y_etet = p_y1[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
        z_etet = p_y1[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

        g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
        g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
        g13 = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;

        c1 = 1/(g11*g33 - g13*g13); // 1/alpha2
        c2 = g13/g33;
        c3 = g13/g11;
        c4 = 1/g22_y1[iptr2];

        temp1 = -c4*((x_xi-c2*x_zt)*x_etet + (y_xi-c2*y_zt)*y_etet + (z_xi-c2*z_zt)*z_etet);
        temp2 = -c4*(x_et*x_etet + y_et*y_etet + z_et*z_etet);
        temp3 = -c4*((x_zt-c3*x_xi)*x_etet + (y_zt-c3*y_xi)*y_etet + (z_zt-c3*z_xi)*z_etet);
        temp_x = g33*x_xixi + g11*x_ztzt - 2*g13*x_xizt;
        temp_y = g33*y_xixi + g11*y_ztzt - 2*g13*y_xizt;
        temp_z = g33*z_xixi + g11*z_ztzt - 2*g13*z_xizt;

        gni = gni1 + (i-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_nx-2) + gni;
        P_y1_loc[iptr] = temp1 - c1*((x_xi-c2*x_zt)*temp_x + (y_xi-c2*y_zt)*temp_y
                       + (z_xi-c2*z_zt)*temp_z);

        Q_y1_loc[iptr] = temp2 - c1*(x_et*temp_x + y_et*temp_y + z_et*temp_z);

        R_y1_loc[iptr] = temp3 - c1*((x_zt-c3*x_xi)*temp_x + (y_zt-c3*y_xi)*temp_y 
                       + (z_zt-c3*z_xi)*temp_z);
      }
    }
  }
  // bdry y2 
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k)
        iptr1 = k*siz_iz + (ny-1)*siz_iy + i+1;   //(i+1,ny-1,k)
        iptr2 = k*siz_iz + (ny-1)*siz_iy + i-1;   //(i-1,ny-1,k)
        iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k+1)
        iptr4 = (k-1)*siz_iz + (ny-1)*siz_iy + i; //(i,ny-1,k-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]); 
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]); 
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]); 
        x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
        y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
        z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
        x_ztzt = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
        y_ztzt = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
        z_ztzt = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

        iptr1 = (k-1)*siz_iz + (ny-1)*siz_iy + i-1; //(i-1,ny-1,k-1)
        iptr2 = (k-1)*siz_iz + (ny-1)*siz_iy + i+1; //(i+1,ny-1,k-1)
        iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i+1; //(i+1,ny-1,k+1)
        iptr4 = (k+1)*siz_iz + (ny-1)*siz_iy + i-1; //(i-1,ny-1,k+1)

        x_xizt = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
        y_xizt = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
        z_xizt = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
        
        iptr1 = k*siz_iz + (ny-2)*siz_iy + i;     //(i,1,k)
        iptr2 = k*nx + i;
        size = nx*nz;
        x_et = x3d[iptr] - x3d[iptr1];
        y_et = y3d[iptr] - y3d[iptr1];
        z_et = z3d[iptr] - z3d[iptr1];
        x_etet = p_y2[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
        y_etet = p_y2[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
        z_etet = p_y2[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

        g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
        g33 = x_zt*x_zt + y_zt*y_zt + z_zt*z_zt;
        g13 = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;

        c1 = 1/(g11*g33 - g13*g13); // 1/alpha2
        c2 = g13/g33;
        c3 = g13/g11;
        c4 = 1/g22_y2[iptr2];

        temp1 = -c4*((x_xi-c2*x_zt)*x_etet + (y_xi-c2*y_zt)*y_etet + (z_xi-c2*z_zt)*z_etet);
        temp2 = -c4*(x_et*x_etet + y_et*y_etet + z_et*z_etet);
        temp3 = -c4*((x_zt-c3*x_xi)*x_etet + (y_zt-c3*y_xi)*y_etet + (z_zt-c3*z_xi)*z_etet);
        temp_x = g33*x_xixi + g11*x_ztzt - 2*g13*x_xizt;
        temp_y = g33*y_xixi + g11*y_ztzt - 2*g13*y_xizt;
        temp_z = g33*z_xixi + g11*z_ztzt - 2*g13*z_xizt;

        gni = gni1 + (i-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_nx-2) + gni;
        P_y2_loc[iptr] = temp1 - c1*((x_xi-c2*x_zt)*temp_x + (y_xi-c2*y_zt)*temp_y
                       + (z_xi-c2*z_zt)*temp_z);

        Q_y2_loc[iptr] = temp2 - c1*(x_et*temp_x + y_et*temp_y + z_et*temp_z);

        R_y2_loc[iptr] = temp3 - c1*((x_zt-c3*x_xi)*temp_x + (y_zt-c3*y_xi)*temp_y 
                       + (z_zt-c3*z_xi)*temp_z);
      }
    }
  }
  // bdry z1 
  if(neighid[4] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = 0*siz_iz + j*siz_iy + i;   //(i,j,0)
        iptr1 = 0*siz_iz + j*siz_iy + i+1;   //(i+1,j,0)
        iptr2 = 0*siz_iz + j*siz_iy + i-1;   //(i-1,j,0)
        iptr3 = 0*siz_iz + (j+1)*siz_iy + i; //(i,j+1,0)
        iptr4 = 0*siz_iz + (j-1)*siz_iy + i; //(i,j-1,0)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
        x_et = 0.5*(x3d[iptr3] - x3d[iptr4]); 
        y_et = 0.5*(y3d[iptr3] - y3d[iptr4]); 
        z_et = 0.5*(z3d[iptr3] - z3d[iptr4]); 
        x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
        y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
        z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
        x_etet = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
        y_etet = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
        z_etet = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

        iptr1 = 0*siz_iz + (j-1)*siz_iy + i-1; //(i-1,j-1,0)
        iptr2 = 0*siz_iz + (j-1)*siz_iy + i+1; //(i+1,j-1,0)
        iptr3 = 0*siz_iz + (j+1)*siz_iy + i+1; //(i+1,j+1,0)
        iptr4 = 0*siz_iz + (j+1)*siz_iy + i-1; //(i-1,j+1,0)

        x_xiet = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
        y_xiet = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
        z_xiet = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
        
        iptr1 = 1*siz_iz + j*siz_iy + i;     //(i,j,1)
        iptr2 = j*nx + i;
        size = nx*ny;
        x_zt = x3d[iptr1] - x3d[iptr];
        y_zt = y3d[iptr1] - y3d[iptr];
        z_zt = z3d[iptr1] - z3d[iptr];
        x_ztzt = p_z1[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
        y_ztzt = p_z1[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
        z_ztzt = p_z1[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

        g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
        g22 = x_et*x_et + y_et*y_et + z_et*z_et;
        g12 = x_xi*x_et + y_xi*y_et + z_xi*z_et;

        c1 = 1/(g11*g22 - g12*g12); // 1/alpha3
        c2 = g12/g22;
        c3 = g12/g11;
        c4 = 1/g33_z1[iptr2];

        temp1 = -c4*((x_xi-c2*x_et)*x_ztzt + (y_xi-c2*y_et)*y_ztzt + (z_xi-c2*z_et)*z_ztzt);
        temp2 = -c4*((x_et-c3*x_xi)*x_ztzt + (y_et-c3*y_xi)*y_ztzt + (z_et-c3*z_xi)*z_ztzt);
        temp3 = -c4*(x_zt*x_ztzt + y_zt*y_ztzt + z_zt*z_ztzt);
        temp_x = g22*x_xixi + g11*x_etet - 2*g12*x_xiet;
        temp_y = g22*y_xixi + g11*y_etet - 2*g12*y_xiet;
        temp_z = g22*z_xixi + g11*z_etet - 2*g12*z_xiet;

        gni = gni1 + (i-1);
        gnj = gnj1 + (j-1);
        iptr = gnj*(total_nx-2) + gni;
        P_z1_loc[iptr] = temp1 - c1*((x_xi-c2*x_et)*temp_x + (y_xi-c2*y_et)*temp_y
                       + (z_xi-c2*z_et)*temp_z);

        Q_z1_loc[iptr] = temp2 - c1*((x_et-c3*x_xi)*temp_x + (y_et-c3*y_xi)*temp_y 
                       + (z_et-c3*z_xi)*temp_z);

        R_z1_loc[iptr] = temp3 - c1*(x_zt*temp_x + y_zt*temp_y + z_zt*temp_z);
      }
    }
  }
  // bdry z2 
  if(neighid[5] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = (nz-1)*siz_iz + j*siz_iy + i;   //(i,j,nz-1)
        iptr1 = (nz-1)*siz_iz + j*siz_iy + i+1;   //(i+1,j,nz-1)
        iptr2 = (nz-1)*siz_iz + j*siz_iy + i-1;   //(i-1,j,nz-1)
        iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i; //(i,j+1,nz-1)
        iptr4 = (nz-1)*siz_iz + (j-1)*siz_iy + i; //(i,j-1,nz-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]); 
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]); 
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]); 
        x_et = 0.5*(x3d[iptr3] - x3d[iptr4]); 
        y_et = 0.5*(y3d[iptr3] - y3d[iptr4]); 
        z_et = 0.5*(z3d[iptr3] - z3d[iptr4]); 
        x_xixi = x3d[iptr1] + x3d[iptr2] - 2*x3d[iptr];
        y_xixi = y3d[iptr1] + y3d[iptr2] - 2*y3d[iptr];
        z_xixi = z3d[iptr1] + z3d[iptr2] - 2*z3d[iptr];
        x_etet = x3d[iptr3] + x3d[iptr4] - 2*x3d[iptr];
        y_etet = y3d[iptr3] + y3d[iptr4] - 2*y3d[iptr];
        z_etet = z3d[iptr3] + z3d[iptr4] - 2*z3d[iptr];

        iptr1 = (nz-1)*siz_iz + (j-1)*siz_iy + i-1; //(i-1,j-1,nz-1)
        iptr2 = (nz-1)*siz_iz + (j-1)*siz_iy + i+1; //(i+1,j-1,nz-1)
        iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i+1; //(i+1,j+1,nz-1)
        iptr4 = (nz-1)*siz_iz + (j+1)*siz_iy + i-1; //(i-1,j+1,nz-1)

        x_xiet = 0.25*(x3d[iptr1] + x3d[iptr3] - x3d[iptr2] - x3d[iptr4]);
        y_xiet = 0.25*(y3d[iptr1] + y3d[iptr3] - y3d[iptr2] - y3d[iptr4]);
        z_xiet = 0.25*(z3d[iptr1] + z3d[iptr3] - z3d[iptr2] - z3d[iptr4]);
        
        iptr1 = (nz-2)*siz_iz + j*siz_iy + i;     //(i,j,nz-2)
        iptr2 = j*nx + i;
        size = nx*ny;
        x_zt = x3d[iptr] - x3d[iptr1];
        y_zt = y3d[iptr] - y3d[iptr1];
        z_zt = z3d[iptr] - z3d[iptr1];
        x_ztzt = p_z2[iptr2+0*size] + x3d[iptr1] - 2*x3d[iptr];
        y_ztzt = p_z2[iptr2+1*size] + y3d[iptr1] - 2*y3d[iptr];
        z_ztzt = p_z2[iptr2+2*size] + z3d[iptr1] - 2*z3d[iptr];

        g11 = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
        g22 = x_et*x_et + y_et*y_et + z_et*z_et;
        g12 = x_xi*x_et + y_xi*y_et + z_xi*z_et;

        c1 = 1/(g11*g22 - g12*g12); // 1/alpha3
        c2 = g12/g22;
        c3 = g12/g11;
        c4 = 1/g33_z2[iptr2];

        temp1 = -c4*((x_xi-c2*x_et)*x_ztzt + (y_xi-c2*y_et)*y_ztzt + (z_xi-c2*z_et)*z_ztzt);
        temp2 = -c4*((x_et-c3*x_xi)*x_ztzt + (y_et-c3*y_xi)*y_ztzt + (z_et-c3*z_xi)*z_ztzt);
        temp3 = -c4*(x_zt*x_ztzt + y_zt*y_ztzt + z_zt*z_ztzt);
        temp_x = g22*x_xixi + g11*x_etet - 2*g12*x_xiet;
        temp_y = g22*y_xixi + g11*y_etet - 2*g12*y_xiet;
        temp_z = g22*z_xixi + g11*z_etet - 2*g12*z_xiet;

        gni = gni1 + (i-1);
        gnj = gnj1 + (j-1);
        iptr = gnj*(total_nx-2) + gni;
        P_z2_loc[iptr] = temp1 - c1*((x_xi-c2*x_et)*temp_x + (y_xi-c2*y_et)*temp_y
                       + (z_xi-c2*z_et)*temp_z);

        Q_z2_loc[iptr] = temp2 - c1*((x_et-c3*x_xi)*temp_x + (y_et-c3*y_xi)*temp_y 
                       + (z_et-c3*z_xi)*temp_z);

        R_z2_loc[iptr] = temp3 - c1*(x_zt*temp_x + y_zt*temp_y + z_zt*temp_z);
      }
    }
  }

  MPI_Barrier(topocomm);

  // not include 2 bdry points
  size_t siz_bx = (total_ny-2)*(total_nz-2);
  size_t siz_by = (total_nx-2)*(total_nz-2);
  size_t siz_bz = (total_nx-2)*(total_ny-2);
  MPI_Allreduce(P_x1_loc, P_x1, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x1_loc, Q_x1, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_x1_loc, R_x1, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_x2_loc, P_x2, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x2_loc, Q_x2, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_x2_loc, R_x2, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(P_y1_loc, P_y1, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_y1_loc, Q_y1, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_y1_loc, R_y1, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_y2_loc, P_y2, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_y2_loc, Q_y2, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_y2_loc, R_y2, siz_by, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(P_z1_loc, P_z1, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z1_loc, Q_z1, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_z1_loc, R_z1, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_z2_loc, P_z2, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z2_loc, Q_z2, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_z2_loc, R_z2, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);

  return 0;
}

int
ghost_cal(float *x3d, float *y3d, float *z3d, int nx,
          int ny, int nz, float *p_x, float *p_y,
          float *p_z, float *g11_x, float *g22_y,
          float *g33_z, int *neighid)
{
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xi0,y_xi0,z_xi0;
  float x_et0,y_et0,z_et0;
  float x_zt0,y_zt0,z_zt0;
  float vn_x,vn_y,vn_z;
  float vn_x0,vn_y0,vn_z0;
  float len_vn,coef;
  float *p_x1,*p_x2;
  float *p_y1,*p_y2;
  float *p_z1,*p_z2;
  float *g11_x1,*g11_x2;
  float *g22_y1,*g22_y2;
  float *g33_z1,*g33_z2;
  size_t siz_iy = nx;
  size_t siz_iz = nx*ny;
  size_t size;
  p_x1 = p_x;
  p_x2 = p_x + ny*nz*3;
  p_y1 = p_y;
  p_y2 = p_y + nx*nz*3;
  p_z1 = p_z;
  p_z2 = p_z + nx*ny*3;
  g11_x1 = g11_x;
  g11_x2 = g11_x + ny*nz;
  g22_y1 = g22_y;
  g22_y2 = g22_y + nx*nz;
  g33_z1 = g33_z;
  g33_z2 = g33_z + nx*ny;
  // bdry x1 r_et X r_zt
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr1 = k*siz_iz + j*siz_iy + 1;    //(1,j,k)
        iptr2 = k*siz_iz + j*siz_iy + 0;    //(0,j,k)
        x_xi0 = x3d[iptr1] - x3d[iptr2];
        y_xi0 = y3d[iptr1] - y3d[iptr2];
        z_xi0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + (j+1)*siz_iy + 0;   //(0,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + 0;   //(0,j-1,k)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr3 = (k+1)*siz_iz + j*siz_iy + 0;   //(0,j,k+1)
        iptr4 = (k-1)*siz_iz + j*siz_iy + 0;   //(0,j,k-1)
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
        // orth vector. cross product
        vn_x = y_et*z_zt - z_et*y_zt;
        vn_y = z_et*x_zt - x_et*z_zt;
        vn_z = x_et*y_zt - y_et*x_zt;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_xi0 to vn. dot product
        coef = x_xi0*vn_x0 + y_xi0*vn_y0 + z_xi0*vn_z0;

        x_xi = coef*vn_x0;
        y_xi = coef*vn_y0;
        z_xi = coef*vn_z0;
         
        iptr = k*ny + j;
        iptr1 = k*siz_iz + j*siz_iy + 0;
        size = ny*nz;
        p_x1[iptr+0*size] = x3d[iptr1] - x_xi;
        p_x1[iptr+1*size] = y3d[iptr1] - y_xi;
        p_x1[iptr+2*size] = z3d[iptr1] - z_xi;
        g11_x1[iptr] = pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2);
      }
    }
  }
  // bdry x2 r_et X r_zt
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr1 = k*siz_iz + j*siz_iy + nx-1;    //(nx-1,j,k)
        iptr2 = k*siz_iz + j*siz_iy + nx-2;    //(nx-2,j,k)
        x_xi0 = x3d[iptr1] - x3d[iptr2];
        y_xi0 = y3d[iptr1] - y3d[iptr2];
        z_xi0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1;   //(nx-1,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1;   //(nx-1,j-1,k)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1;   //(nx-1,j,k+1)
        iptr4 = (k-1)*siz_iz + j*siz_iy + nx-1;   //(nx-1,j,k-1)
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
        // orth vector. cross product
        vn_x = y_et*z_zt - z_et*y_zt;
        vn_y = z_et*x_zt - x_et*z_zt;
        vn_z = x_et*y_zt - y_et*x_zt;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_xi0 to vn. dot product
        coef = x_xi0*vn_x0 + y_xi0*vn_y0 + z_xi0*vn_z0;

        x_xi = coef*vn_x0;
        y_xi = coef*vn_y0;
        z_xi = coef*vn_z0;
         
        iptr = k*ny + j;
        iptr1 = k*siz_iz + j*siz_iy + nx-1;
        size = ny*nz;
        p_x2[iptr+0*size] = x3d[iptr1] + x_xi;
        p_x2[iptr+1*size] = y3d[iptr1] + y_xi;
        p_x2[iptr+2*size] = z3d[iptr1] + z_xi;
        g11_x2[iptr] = pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2);
      }
    }
  }
  // bdry y1 r_zt X r_xi
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = k*siz_iz + 1*siz_iy + i;    //(i,1,k)
        iptr2 = k*siz_iz + 0*siz_iy + i;    //(i,0,k)
        x_et0 = x3d[iptr1] - x3d[iptr2];
        y_et0 = y3d[iptr1] - y3d[iptr2];
        z_et0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + 0*siz_iy + (i+1);   //(i+1,0,k)
        iptr2 = k*siz_iz + 0*siz_iy + (i-1);   //(i-1,0,k)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr3 = (k+1)*siz_iz + 0*siz_iy + i;   //(i,0,k+1)
        iptr4 = (k-1)*siz_iz + 0*siz_iy + i;   //(i,0,k-1)
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
        // orth vector. cross product
        vn_x = y_zt*z_xi - z_zt*y_xi;
        vn_y = z_zt*x_xi - x_zt*z_xi;
        vn_z = x_zt*y_xi - y_zt*x_xi;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        coef = x_et0*vn_x0 + y_et0*vn_y0 + z_et0*vn_z0;

        x_et = coef*vn_x0;
        y_et = coef*vn_y0;
        z_et = coef*vn_z0;
         
        iptr = k*nx + i;
        iptr1 = k*siz_iz + 0*siz_iy + i;
        size =  nx*nz;
        p_y1[iptr+0*size] = x3d[iptr1] - x_et;
        p_y1[iptr+1*size] = y3d[iptr1] - y_et;
        p_y1[iptr+2*size] = z3d[iptr1] - z_et;
        g22_y1[iptr] = pow(x_et,2) + pow(y_et,2) + pow(z_et,2);
      }
    }
  }
  // bdry y2 r_zt X r_xi
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = k*siz_iz + (ny-1)*siz_iy + i;    //(i,ny-1,k)
        iptr2 = k*siz_iz + (ny-2)*siz_iy + i;    //(i,ny-2,k)
        x_et0 = x3d[iptr1] - x3d[iptr2];
        y_et0 = y3d[iptr1] - y3d[iptr2];
        z_et0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + (ny-1)*siz_iy + (i+1);   //(i+1,ny-1,k)
        iptr2 = k*siz_iz + (ny-1)*siz_iy + (i-1);   //(i-1,ny-1,k)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k+1)
        iptr4 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k-1)
        x_zt = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3] - z3d[iptr4]);
        // orth vector. cross product
        vn_x = y_zt*z_xi - z_zt*y_xi;
        vn_y = z_zt*x_xi - x_zt*z_xi;
        vn_z = x_zt*y_xi - y_zt*x_xi;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        coef = x_et0*vn_x0 + y_et0*vn_y0 + z_et0*vn_z0;

        x_et = coef*vn_x0;
        y_et = coef*vn_y0;
        z_et = coef*vn_z0;
         
        iptr = k*nx + i;
        iptr1 = k*siz_iz + (ny-1)*siz_iy + i;
        size = nx*nz;
        p_y2[iptr+0*size] = x3d[iptr1] + x_et;
        p_y2[iptr+1*size] = y3d[iptr1] + y_et;
        p_y2[iptr+2*size] = z3d[iptr1] + z_et;
        g22_y2[iptr] = pow(x_et,2) + pow(y_et,2) + pow(z_et,2);
      }
    }
  }
  // bdry z1 r_xi X r_et
  if(neighid[4] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = 1*siz_iz + j*siz_iy + i;    //(i,j,1)
        iptr2 = 0*siz_iz + j*siz_iy + i;    //(i,j,0)
        x_zt0 = x3d[iptr1] - x3d[iptr2];
        y_zt0 = y3d[iptr1] - y3d[iptr2];
        z_zt0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = 0*siz_iz + j*siz_iy + (i+1);   //(i+1,j,0)
        iptr2 = 0*siz_iz + j*siz_iy + (i-1);   //(i-1,j,0)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr3 = 0*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,0)
        iptr4 = 0*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,0)
        x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);
        // orth vector. cross product
        vn_x = y_xi*z_et - z_xi*y_et;
        vn_y = z_xi*x_et - x_xi*z_et;
        vn_z = x_xi*y_et - y_xi*x_et;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        coef = x_zt0*vn_x0 + y_zt0*vn_y0 + z_zt0*vn_z0;

        x_zt = coef*vn_x0;
        y_zt = coef*vn_y0;
        z_zt = coef*vn_z0;
         
        iptr = j*nx + i;
        iptr1 = 0*siz_iz + j*siz_iy + i;
        size = nx*ny;
        p_z1[iptr+0*size] = x3d[iptr1] - x_zt;
        p_z1[iptr+1*size] = y3d[iptr1] - y_zt;
        p_z1[iptr+2*size] = z3d[iptr1] - z_zt;
        g33_z1[iptr] = pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2);
      }
    }
  }
  // bdry z2 r_xi X r_et
  if(neighid[5] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = (nz-1)*siz_iz + j*siz_iy + i;    //(i,j,nz-1)
        iptr2 = (nz-2)*siz_iz + j*siz_iy + i;    //(i,j,nz-2)
        x_zt0 = x3d[iptr1] - x3d[iptr2];
        y_zt0 = y3d[iptr1] - y3d[iptr2];
        z_zt0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = (nz-1)*siz_iz + j*siz_iy + (i+1);   //(i+1,j,nz-1)
        iptr2 = (nz-1)*siz_iz + j*siz_iy + (i-1);   //(i-1,j,nz-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,nz-1)
        iptr4 = (nz-1)*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,nz-1)
        x_et = 0.5*(x3d[iptr3] - x3d[iptr4]);
        y_et = 0.5*(y3d[iptr3] - y3d[iptr4]);
        z_et = 0.5*(z3d[iptr3] - z3d[iptr4]);
        // orth vector. cross product
        vn_x = y_xi*z_et - z_xi*y_et;
        vn_y = z_xi*x_et - x_xi*z_et;
        vn_z = x_xi*y_et - y_xi*x_et;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        coef = x_zt0*vn_x0 + y_zt0*vn_y0 + z_zt0*vn_z0;

        x_zt = coef*vn_x0;
        y_zt = coef*vn_y0;
        z_zt = coef*vn_z0;
         
        iptr = j*nx + i;
        iptr1 = (nz-1)*siz_iz + j*siz_iy + i;
        size = nx*ny;
        p_z2[iptr+0*size] = x3d[iptr1] + x_zt;
        p_z2[iptr+1*size] = y3d[iptr1] + y_zt;
        p_z2[iptr+2*size] = z3d[iptr1] + z_zt;
        g33_z2[iptr] = pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2);
      }
    }
  }

  return 0;
}

int
higen_gene(gd_t *gdcurv, par_t *par, mympi_t *mympi)
{
  float err = par->iter_err;
  int max_iter = par->max_iter;
  float *coef = par->coef;

  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  size_t iptr, iptr1, iptr2, iptr3;

  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  float *x3d_tmp = gdcurv->x3d_tmp;
  float *y3d_tmp = gdcurv->y3d_tmp;
  float *z3d_tmp = gdcurv->z3d_tmp;

  MPI_Comm comm = mympi->comm;
  int myid = mympi->myid;
  int *neighid = mympi->neighid;

  int num_of_s_reqs = 6;
  int num_of_r_reqs = 6;
  
  src_t *src = (src_t *) malloc(sizeof(src_t)); 
  init_src(src,gdcurv);

  float *dx1 = (float *)mem_calloc_1d_float(nz*ny, 0.0, "dx1");
  float *dx2 = (float *)mem_calloc_1d_float(nz*ny, 0.0, "dx2");
  float *dy1 = (float *)mem_calloc_1d_float(nz*nx, 0.0, "dy1");
  float *dy2 = (float *)mem_calloc_1d_float(nz*nx, 0.0, "dy2");
  float *dz1 = (float *)mem_calloc_1d_float(ny*nx, 0.0, "dz1");
  float *dz2 = (float *)mem_calloc_1d_float(ny*nx, 0.0, "dz2");
  dist_cal(gdcurv,dx1,dx2,dy1,dy2,dz1,dz2,neighid);

  set_src_higen(x3d,y3d,z3d,gdcurv,src,dx1,dx2,dy1,dy2,dz1,dz2,mympi);
  interp_inner_source(src, gdcurv, coef);

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*siz_iz + j*siz_iy + i;
        x3d_tmp[iptr] = x3d[iptr];
        y3d_tmp[iptr] = y3d[iptr];
        z3d_tmp[iptr] = z3d[iptr];
      }
    }
  }

  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  float *ptr_x, *ptr_y, *ptr_z;
  
  int flag_true = 1;
  float resi, resj, resk;
  float max_resi, max_resj, max_resk;
  float dif, dif1, dif2, dif3, dif_x, dif_y, dif_z;
  while(flag_true)
  {
    // update solver
    update_SOR(x3d,y3d,z3d,x3d_tmp,y3d_tmp,z3d_tmp,nx,ny,nz,src->P,src->Q,src->R,w);
    Niter += 1;

    MPI_Startall(num_of_r_reqs, mympi->r_reqs);

    grid_pack_mesg(mympi,gdcurv,x3d_tmp,y3d_tmp,z3d_tmp);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs);

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs,MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs,MPI_STATUS_IGNORE);

    grid_unpack_mesg(mympi,gdcurv,x3d_tmp,y3d_tmp,z3d_tmp);

    max_resi = 0.0;
    max_resj = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++) {
        for(int i=1; i<nx-1; i++)
        {
          iptr  = k*siz_iz + j*siz_iy + i;
          iptr1 = k*siz_iz + j*siz_iy + i+1;
          iptr2 = k*siz_iz + (j+1)*siz_iy  + i;
          iptr3 = (k+1)*siz_iz + j*siz_iy  + i;
          dif_x = x3d_tmp[iptr]-x3d[iptr];
          dif_y = y3d_tmp[iptr]-y3d[iptr];
          dif_z = z3d_tmp[iptr]-z3d[iptr];
          dif = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          dif_x = x3d[iptr1]-x3d[iptr];
          dif_y = y3d[iptr1]-y3d[iptr];
          dif_z = z3d[iptr1]-z3d[iptr];
          dif1 = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          dif_x = x3d[iptr2]-x3d[iptr];
          dif_y = y3d[iptr2]-y3d[iptr];
          dif_z = z3d[iptr2]-z3d[iptr];
          dif2 = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          dif_x = x3d[iptr3]-x3d[iptr];
          dif_y = y3d[iptr3]-y3d[iptr];
          dif_z = z3d[iptr3]-z3d[iptr];
          dif3 = sqrt(pow(dif_x,2) + pow(dif_y,2) + pow(dif_z,2));
          resi = dif/dif1;
          resj = dif/dif2;
          resk = dif/dif3;
          max_resi = fmax(max_resi,resi);
          max_resj = fmax(max_resj,resj);
          max_resk = fmax(max_resk,resk);
        }
      }
    }

    float sendbuf = max_resi;
    MPI_Allreduce(&sendbuf, &max_resi, 1, MPI_REAL, MPI_MAX, comm);
    sendbuf = max_resj;
    MPI_Allreduce(&sendbuf, &max_resj, 1, MPI_REAL, MPI_MAX, comm);
    sendbuf = max_resk;
    MPI_Allreduce(&sendbuf, &max_resk, 1, MPI_REAL, MPI_MAX, comm);

    if(myid == 0)
    {
      fprintf(stdout,"number of iter is %d\n", Niter);
      fprintf(stdout,"max_resi is %f, max_resj is %f, max_resk is %f\n", max_resi, max_resj, max_resk);
    }

    // swap pointer, avoid coping
    ptr_x = x3d;
    ptr_y = y3d;
    ptr_z = z3d;
    x3d = x3d_tmp;
    y3d = y3d_tmp;
    z3d = z3d_tmp;
    x3d_tmp = ptr_x;  
    y3d_tmp = ptr_y;
    z3d_tmp = ptr_z;
    
    set_src_higen(x3d,y3d,z3d,gdcurv,src,dx1,dx2,dy1,dy2,dz1,dz2,mympi);
    interp_inner_source(src, gdcurv, coef);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resj < err && max_resk < err)
    {
      flag_true = 0;
    }

  }

  gdcurv->v4d = x3d;
  gdcurv->x3d = x3d;
  gdcurv->y3d = y3d;
  gdcurv->z3d = z3d;

  gdcurv->v4d_tmp = x3d_tmp;

  free(gdcurv->v4d_tmp);  // free temp grid space
  free(dx1);  
  free(dx2);  
  free(dy1);  
  free(dy2);  
  free(dz1);  
  free(dz2);  

  return 0;
}

int
set_src_higen(float *x3d, float *y3d, float *z3d, 
              gd_t *gdcurv, src_t *src, float *dx1, 
              float *dx2, float *dy1, float *dy2,
              float *dz1, float *dz2, mympi_t *mympi)
             
{
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int ni = gdcurv->ni;
  int nj = gdcurv->nj;
  int nk = gdcurv->nk;
  int gni1 = gdcurv->gni1;
  int gnj1 = gdcurv->gnj1;
  int gnk1 = gdcurv->gnk1;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  int gni,gnj,gnk;

  MPI_Comm topocomm = mympi->topocomm;
  int *neighid = mympi->neighid;

  float *P_x1_loc = src->P_x1_loc;
  float *Q_x1_loc = src->Q_x1_loc;
  float *R_x1_loc = src->R_x1_loc;
  float *P_x2_loc = src->P_x2_loc;
  float *Q_x2_loc = src->Q_x2_loc;
  float *R_x2_loc = src->R_x2_loc;

  float *P_y1_loc = src->P_y1_loc;
  float *Q_y1_loc = src->Q_y1_loc;
  float *R_y1_loc = src->R_y1_loc;
  float *P_y2_loc = src->P_y2_loc;
  float *Q_y2_loc = src->Q_y2_loc;
  float *R_y2_loc = src->R_y2_loc;
  
  float *P_z1_loc = src->P_z1_loc;
  float *Q_z1_loc = src->Q_z1_loc;
  float *R_z1_loc = src->R_z1_loc;
  float *P_z2_loc = src->P_z2_loc;
  float *Q_z2_loc = src->Q_z2_loc;
  float *R_z2_loc = src->R_z2_loc;

  float *P_x1 = src->P_x1;
  float *Q_x1 = src->Q_x1;
  float *R_x1 = src->R_x1;
  float *P_x2 = src->P_x2;
  float *Q_x2 = src->Q_x2;
  float *R_x2 = src->R_x2;

  float *P_y1 = src->P_y1;
  float *Q_y1 = src->Q_y1;
  float *R_y1 = src->R_y1;
  float *P_y2 = src->P_y2;
  float *Q_y2 = src->Q_y2;
  float *R_y2 = src->R_y2;
  
  float *P_z1 = src->P_z1;
  float *Q_z1 = src->Q_z1;
  float *R_z1 = src->R_z1;
  float *P_z2 = src->P_z2;
  float *Q_z2 = src->Q_z2;
  float *R_z2 = src->R_z2;

  float theta0 = PI/2;
  // test a value, not effct iter number
  float a = 0.2; 
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float dot, len_xi, len_et, len_zt, dif_dis;
  float cos_theta, theta, dif_theta;
  float x_xi, y_xi, z_xi;
  float x_et, y_et, z_et;
  float x_zt, y_zt, z_zt;

  // NOTE: P Q R update need change sign
  // when xi = 1, et = 1, zt = 1. 
  // for example:
  // xi = 0. Q R -, P + 
  // xi = 1. Q R +, P - 
  // bdry x1 xi=0
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr  = k*siz_iz + j*siz_iy + 0;         //(0,j,k)
        iptr1 = k*siz_iz + (j+1)*siz_iy + 0;     //(0,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + 0;     //(0,j-1,k)
        iptr3 = (k+1)*siz_iz + j*siz_iy + 0;     //(0,j,k+1)
        iptr4 = (k-1)*siz_iz + j*siz_iy + 0;     //(0,j,k-1)

        x_et = 0.5*(x3d[iptr1]-x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1]-y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1]-z3d[iptr2]);
        x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

        iptr1 = k*siz_iz + j*siz_iy + 1;     //(1,j,k)

        x_xi = x3d[iptr1] - x3d[iptr];
        y_xi = y3d[iptr1] - y3d[iptr];
        z_xi = z3d[iptr1] - z3d[iptr];
        
        gnj = gnj1 + (j-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_ny-2) + gnj;

        // cos(theta) = a.b/(|a|*|b|)
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        // r_xi . r_et
        dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
        cos_theta = dot/(len_xi*len_et);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        Q_x1_loc[iptr] = Q_x1_loc[iptr] - a*tanh(dif_theta);

        // r_xi . r_zt
        dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
        cos_theta = dot/(len_xi*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        R_x1_loc[iptr] = R_x1_loc[iptr] - a*tanh(dif_theta);

        iptr1 = k*ny+j;
        dif_dis = (dx1[iptr1]-len_xi)/dx1[iptr1];
        P_x1_loc[iptr] = P_x1_loc[iptr] + a*tanh(dif_dis);
      }
    }
  }
  // bdry x2 xi=1
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr  = k*siz_iz + j*siz_iy + nx-1;         //(nx-1,j,k)
        iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1;     //(nx-1,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1;     //(nx-1,j-1,k)
        iptr3 = (k+1)*siz_iz + j*siz_iy + nx-1;     //(nx-1,j,k+1)
        iptr4 = (k-1)*siz_iz + j*siz_iy + nx-1;     //(nx-1,j,k-1)

        x_et = 0.5*(x3d[iptr1]-x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1]-y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1]-z3d[iptr2]);
        x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

        iptr1 = k*siz_iz + j*siz_iy + nx-2;     //(nx-2,j,k)

        x_xi = x3d[iptr] - x3d[iptr1];
        y_xi = y3d[iptr] - y3d[iptr1];
        z_xi = z3d[iptr] - z3d[iptr1];
        
        gnj = gnj1 + (j-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_ny-2) + gnj;

        // cos(theta) = a.b/(|a|*|b|)
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        // r_xi . r_et
        dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
        cos_theta = dot/(len_xi*len_et);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        Q_x2_loc[iptr] = Q_x2_loc[iptr] + a*tanh(dif_theta);

        // r_xi . r_zt
        dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
        cos_theta = dot/(len_xi*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        R_x2_loc[iptr] = R_x2_loc[iptr] + a*tanh(dif_theta);

        iptr1 = k*ny+j;
        dif_dis = (dx2[iptr1]-len_xi)/dx2[iptr1];
        P_x2_loc[iptr] = P_x2_loc[iptr] - a*tanh(dif_dis);
      }
    }
  }
  // bdry y1 et=0
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*siz_iz + 0*siz_iy + i;       //(i,0,k)
        iptr1 = k*siz_iz + 0*siz_iy + i+1;     //(i+1,0,k)
        iptr2 = k*siz_iz + 0*siz_iy + i-1;     //(i-1,0,k)
        iptr3 = (k+1)*siz_iz + 0*siz_iy + i;   //(i,0,k+1)
        iptr4 = (k-1)*siz_iz + 0*siz_iy + i;   //(i,0,k-1)

        x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
        x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

        iptr1 = k*siz_iz + 1*siz_iy + i;     //(i,1,k)

        x_et = x3d[iptr1] - x3d[iptr];
        y_et = y3d[iptr1] - y3d[iptr];
        z_et = z3d[iptr1] - z3d[iptr];
        
        gni = gni1 + (i-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_nx-2) + gni;

        // cos(theta) = a.b/(|a|*|b|)
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        // r_xi . r_et
        dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
        cos_theta = dot/(len_xi*len_et);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        P_y1_loc[iptr] = P_y1_loc[iptr] - a*tanh(dif_theta);

        // r_et . r_zt
        dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
        cos_theta = dot/(len_et*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        R_y1_loc[iptr] = R_y1_loc[iptr] - a*tanh(dif_theta);

        iptr1 = k*nx+i;
        dif_dis = (dy1[iptr1]-len_et)/dy1[iptr1];
        Q_y1_loc[iptr] = Q_y1_loc[iptr] + a*tanh(dif_dis);
      }
    }
  }
  // bdry y2 et=0
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*siz_iz + (ny-1)*siz_iy + i;       //(i,ny-1,k)
        iptr1 = k*siz_iz + (ny-1)*siz_iy + i+1;     //(i+1,ny-1,k)
        iptr2 = k*siz_iz + (ny-1)*siz_iy + i-1;     //(i-1,ny-1,k)
        iptr3 = (k+1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k+1)
        iptr4 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k-1)

        x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
        x_zt = 0.5*(x3d[iptr3]-x3d[iptr4]);
        y_zt = 0.5*(y3d[iptr3]-y3d[iptr4]);
        z_zt = 0.5*(z3d[iptr3]-z3d[iptr4]);

        iptr1 = k*siz_iz + (ny-2)*siz_iy + i;     //(i,ny-2,k)

        x_et = x3d[iptr] - x3d[iptr1];
        y_et = y3d[iptr] - y3d[iptr1];
        z_et = z3d[iptr] - z3d[iptr1];
        
        gni = gni1 + (i-1);
        gnk = gnk1 + (k-1);
        iptr = gnk*(total_nx-2) + gni;

        // cos(theta) = a.b/(|a|*|b|)
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        // r_xi . r_et
        dot = x_xi*x_et + y_xi*y_et + z_xi*z_et;
        cos_theta = dot/(len_xi*len_et);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        P_y2_loc[iptr] = P_y2_loc[iptr] + a*tanh(dif_theta);

        // r_et . r_zt
        dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
        cos_theta = dot/(len_et*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        R_y2_loc[iptr] = R_y2_loc[iptr] + a*tanh(dif_theta);

        iptr1 = k*nx+i;
        dif_dis = (dy2[iptr1]-len_et)/dy2[iptr1];
        Q_y2_loc[iptr] = Q_y2_loc[iptr] - a*tanh(dif_dis);
      }
    }
  }
  // bdry z1 zt=0
  if(neighid[4] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = 0*siz_iz + j*siz_iy + i;       //(i,j,0)
        iptr1 = 0*siz_iz + j*siz_iy + i+1;     //(i+1,j,0)
        iptr2 = 0*siz_iz + j*siz_iy + i-1;     //(i-1,j,0)
        iptr3 = 0*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,0)
        iptr4 = 0*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,0)

        x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
        x_et = 0.5*(x3d[iptr3]-x3d[iptr4]);
        y_et = 0.5*(y3d[iptr3]-y3d[iptr4]);
        z_et = 0.5*(z3d[iptr3]-z3d[iptr4]);

        iptr1 = 1*siz_iz + j*siz_iy + i;     //(i,j,1)

        x_zt = x3d[iptr1] - x3d[iptr];
        y_zt = y3d[iptr1] - y3d[iptr];
        z_zt = z3d[iptr1] - z3d[iptr];
        
        gni = gni1 + (i-1);
        gnj = gnj1 + (j-1);
        iptr = gnj*(total_nx-2) + gni;

        // cos(theta) = a.b/(|a|*|b|)
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        // r_xi . r_zt
        dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
        cos_theta = dot/(len_xi*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        P_z1_loc[iptr] = P_z1_loc[iptr] - a*tanh(dif_theta);

        // r_et . r_zt
        dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
        cos_theta = dot/(len_et*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        Q_z1_loc[iptr] = Q_z1_loc[iptr] - a*tanh(dif_theta);

        iptr1 = j*nx+i;
        dif_dis = (dz1[iptr1]-len_zt)/dz1[iptr1];
        R_z1_loc[iptr] = R_z1_loc[iptr] + a*tanh(dif_dis);
      }
    }
  }
  // bdry z2 zt=1
  if(neighid[5] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = (nz-1)*siz_iz + j*siz_iy + i;       //(i,j,nz-1)
        iptr1 = (nz-1)*siz_iz + j*siz_iy + i+1;     //(i+1,j,nz-1)
        iptr2 = (nz-1)*siz_iz + j*siz_iy + i-1;     //(i-1,j,nz-1)
        iptr3 = (nz-1)*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,nz-1)
        iptr4 = (nz-1)*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,nz-1)

        x_xi = 0.5*(x3d[iptr1]-x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1]-y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1]-z3d[iptr2]);
        x_et = 0.5*(x3d[iptr3]-x3d[iptr4]);
        y_et = 0.5*(y3d[iptr3]-y3d[iptr4]);
        z_et = 0.5*(z3d[iptr3]-z3d[iptr4]);

        iptr1 = (nz-2)*siz_iz + j*siz_iy + i;     //(i,j,nz-2)

        x_zt = x3d[iptr] - x3d[iptr1];
        y_zt = y3d[iptr] - y3d[iptr1];
        z_zt = z3d[iptr] - z3d[iptr1];
        
        gni = gni1 + (i-1);
        gnj = gnj1 + (j-1);
        iptr = gnj*(total_nx-2) + gni;

        // cos(theta) = a.b/(|a|*|b|)
        len_xi = sqrt(pow(x_xi,2) + pow(y_xi,2) + pow(z_xi,2));
        len_et = sqrt(pow(x_et,2) + pow(y_et,2) + pow(z_et,2));
        len_zt = sqrt(pow(x_zt,2) + pow(y_zt,2) + pow(z_zt,2));
        // r_xi . r_zt
        dot = x_xi*x_zt + y_xi*y_zt + z_xi*z_zt;
        cos_theta = dot/(len_xi*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        P_z2_loc[iptr] = P_z2_loc[iptr] + a*tanh(dif_theta);

        // r_et . r_zt
        dot = x_et*x_zt + y_et*y_zt + z_et*z_zt;
        cos_theta = dot/(len_et*len_zt);
        theta = acos(cos_theta);
        dif_theta = (theta0-theta)/theta0;
        Q_z2_loc[iptr] = Q_z2_loc[iptr] + a*tanh(dif_theta);

        iptr1 = j*nx+i;
        dif_dis = (dz2[iptr1]-len_zt)/dz2[iptr1];
        R_z2_loc[iptr] = R_z2_loc[iptr] - a*tanh(dif_dis);
      }
    }
  }

  MPI_Barrier(topocomm);

  // not include 2 bdry points
  size_t siz_bx = (total_ny-2)*(total_nz-2);
  size_t siz_by = (total_nx-2)*(total_nz-2);
  size_t siz_bz = (total_nx-2)*(total_ny-2);
  MPI_Allreduce(P_x1_loc, P_x1, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x1_loc, Q_x1, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_x1_loc, R_x1, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_x2_loc, P_x2, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x2_loc, Q_x2, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_x2_loc, R_x2, siz_bx, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(P_y1_loc, P_y1, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_y1_loc, Q_y1, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_y1_loc, R_y1, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_y2_loc, P_y2, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_y2_loc, Q_y2, siz_by, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_y2_loc, R_y2, siz_by, MPI_FLOAT, MPI_SUM, topocomm);

  MPI_Allreduce(P_z1_loc, P_z1, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z1_loc, Q_z1, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_z1_loc, R_z1, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_z2_loc, P_z2, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z2_loc, Q_z2, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(R_z2_loc, R_z2, siz_bz, MPI_FLOAT, MPI_SUM, topocomm);

  return 0;
}

int
dist_cal(gd_t *gdcurv, float *dx1, float *dx2, float *dy1, float *dy2, float *dz1, float *dz2, int *neighid)
{
  size_t iptr,iptr1,iptr2;
  float x_xi,y_xi,z_xi;
  float x_et,y_et,z_et;
  float x_zt,y_zt,z_zt;
  float x_xi0,y_xi0,z_xi0;
  float x_et0,y_et0,z_et0;
  float x_zt0,y_zt0,z_zt0;
  float vn_x,vn_y,vn_z,len_vn;
  float vn_x0,vn_y0,vn_z0;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;
  float *x3d = gdcurv->x3d; 
  float *y3d = gdcurv->y3d; 
  float *z3d = gdcurv->z3d; 
  // bdry x1 r_et X r_zt
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr1 = k*siz_iz + j*siz_iy + 1;    //(1,j,k)
        iptr2 = k*siz_iz + j*siz_iy + 0;    //(0,j,k)
        x_xi0 = x3d[iptr1] - x3d[iptr2];
        y_xi0 = y3d[iptr1] - y3d[iptr2];
        z_xi0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + (j+1)*siz_iy + 0;   //(0,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + 0;   //(0,j-1,k)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr1 = (k+1)*siz_iz + j*siz_iy + 0;   //(0,j,k+1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + 0;   //(0,j,k-1)
        x_zt = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_zt = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_zt = 0.5*(z3d[iptr1] - z3d[iptr2]);
        // orth vector. cross product
        vn_x = y_et*z_zt - z_et*y_zt;
        vn_y = z_et*x_zt - x_et*z_zt;
        vn_z = x_et*y_zt - y_et*x_zt;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_xi0 to vn. dot product
        iptr = k*ny + j;
        dx1[iptr] = x_xi0*vn_x0 + y_xi0*vn_y0 + z_xi0*vn_z0;
      }
    }
  }
  // bdry x2 r_et X r_zt
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int j=1; j<ny-1; j++)
      {
        iptr1 = k*siz_iz + j*siz_iy + nx-1;    //(nx-1,j,k)
        iptr2 = k*siz_iz + j*siz_iy + nx-2;    //(nx-2,j,k)
        x_xi0 = x3d[iptr1] - x3d[iptr2];
        y_xi0 = y3d[iptr1] - y3d[iptr2];
        z_xi0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + (j+1)*siz_iy + nx-1;   //(nx-1,j+1,k)
        iptr2 = k*siz_iz + (j-1)*siz_iy + nx-1;   //(nx-1,j-1,k)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr1 = (k+1)*siz_iz + j*siz_iy + nx-1;   //(nx-1,j,k+1)
        iptr2 = (k-1)*siz_iz + j*siz_iy + nx-1;   //(nx-1,j,k-1)
        x_zt = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_zt = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_zt = 0.5*(z3d[iptr1] - z3d[iptr2]);
        // orth vector. cross product
        vn_x = y_et*z_zt - z_et*y_zt;
        vn_y = z_et*x_zt - x_et*z_zt;
        vn_z = x_et*y_zt - y_et*x_zt;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_xi0 to vn. dot product
        iptr = k*ny + j;
        dx2[iptr] = x_xi0*vn_x0 + y_xi0*vn_y0 + z_xi0*vn_z0;
      }
    }
  }
  // bdry y1 r_zt X r_xi
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = k*siz_iz + 1*siz_iy + i;    //(i,1,k)
        iptr2 = k*siz_iz + 0*siz_iy + i;    //(i,0,k)
        x_et0 = x3d[iptr1] - x3d[iptr2];
        y_et0 = y3d[iptr1] - y3d[iptr2];
        z_et0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + 0*siz_iy + (i+1);   //(i+1,0,k)
        iptr2 = k*siz_iz + 0*siz_iy + (i-1);   //(i-1,0,k)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr1 = (k+1)*siz_iz + 0*siz_iy + i;   //(i,0,k+1)
        iptr2 = (k-1)*siz_iz + 0*siz_iy + i;   //(i,0,k-1)
        x_zt = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_zt = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_zt = 0.5*(z3d[iptr1] - z3d[iptr2]);
        // orth vector. cross product
        vn_x = y_zt*z_xi - z_zt*y_xi;
        vn_y = z_zt*x_xi - x_zt*z_xi;
        vn_z = x_zt*y_xi - y_zt*x_xi;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        iptr = k*nx + i;
        dy1[iptr] = x_et0*vn_x0 + y_et0*vn_y0 + z_et0*vn_z0;
      }
    }
  }
  // bdry y2 r_zt X r_xi
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = k*siz_iz + (ny-1)*siz_iy + i;    //(i,ny-1,k)
        iptr2 = k*siz_iz + (ny-2)*siz_iy + i;    //(i,ny-2,k)
        x_et0 = x3d[iptr1] - x3d[iptr2];
        y_et0 = y3d[iptr1] - y3d[iptr2];
        z_et0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = k*siz_iz + (ny-1)*siz_iy + (i+1);   //(i+1,ny-1,k)
        iptr2 = k*siz_iz + (ny-1)*siz_iy + (i-1);   //(i-1,ny-1,k)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr1 = (k+1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k+1)
        iptr2 = (k-1)*siz_iz + (ny-1)*siz_iy + i;   //(i,ny-1,k-1)
        x_zt = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_zt = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_zt = 0.5*(z3d[iptr1] - z3d[iptr2]);
        // orth vector. cross product
        vn_x = y_zt*z_xi - z_zt*y_xi;
        vn_y = z_zt*x_xi - x_zt*z_xi;
        vn_z = x_zt*y_xi - y_zt*x_xi;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        iptr = k*nx + i;
        dy2[iptr] = x_et0*vn_x0 + y_et0*vn_y0 + z_et0*vn_z0;
      }
    }
  }
  // bdry z1 r_xi X r_et
  if(neighid[4] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = 1*siz_iz + j*siz_iy + i;    //(i,j,1)
        iptr2 = 0*siz_iz + j*siz_iy + i;    //(i,j,0)
        x_zt0 = x3d[iptr1] - x3d[iptr2];
        y_zt0 = y3d[iptr1] - y3d[iptr2];
        z_zt0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = 0*siz_iz + j*siz_iy + (i+1);   //(i+1,j,0)
        iptr2 = 0*siz_iz + j*siz_iy + (i-1);   //(i-1,j,0)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr1 = 0*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,0)
        iptr2 = 0*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,0)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        // orth vector. cross product
        vn_x = y_xi*z_et - z_xi*y_et;
        vn_y = z_xi*x_et - x_xi*z_et;
        vn_z = x_xi*y_et - y_xi*x_et;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        iptr = j*nx + i;
        dz1[iptr] = x_zt0*vn_x0 + y_zt0*vn_y0 + z_zt0*vn_z0;
      }
    }
  }
  // bdry z2 r_xi X r_et
  if(neighid[5] == MPI_PROC_NULL)
  {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr1 = (nz-1)*siz_iz + j*siz_iy + i;    //(i,j,nz-1)
        iptr2 = (nz-2)*siz_iz + j*siz_iy + i;    //(i,j,nz-2)
        x_zt0 = x3d[iptr1] - x3d[iptr2];
        y_zt0 = y3d[iptr1] - y3d[iptr2];
        z_zt0 = z3d[iptr1] - z3d[iptr2];
        iptr1 = (nz-1)*siz_iz + j*siz_iy + (i+1);   //(i+1,j,nz-1)
        iptr2 = (nz-1)*siz_iz + j*siz_iy + (i-1);   //(i-1,j,nz-1)
        x_xi = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_xi = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_xi = 0.5*(z3d[iptr1] - z3d[iptr2]);
        iptr1 = (nz-1)*siz_iz + (j+1)*siz_iy + i;   //(i,j+1,nz-1)
        iptr2 = (nz-1)*siz_iz + (j-1)*siz_iy + i;   //(i,j-1,nz-1)
        x_et = 0.5*(x3d[iptr1] - x3d[iptr2]);
        y_et = 0.5*(y3d[iptr1] - y3d[iptr2]);
        z_et = 0.5*(z3d[iptr1] - z3d[iptr2]);
        // orth vector. cross product
        vn_x = y_xi*z_et - z_xi*y_et;
        vn_y = z_xi*x_et - x_xi*z_et;
        vn_z = x_xi*y_et - y_xi*x_et;
        len_vn = sqrt(pow(vn_x,2)+pow(vn_y,2)+pow(vn_z,2));
        // norm
        vn_x0 = vn_x/len_vn;
        vn_y0 = vn_y/len_vn;
        vn_z0 = vn_z/len_vn;
        // projection from r_et0 to vn. dot product
        iptr = j*nx + i;
        dz2[iptr] = x_zt0*vn_x0 + y_zt0*vn_y0 + z_zt0*vn_z0;
      }
    }
  }

  return 0;
}

int interp_inner_source(src_t *src, gd_t *gdcurv, float *coef)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nj1 = gdcurv->nj1;
  int nj2 = gdcurv->nj2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx = gdcurv->nx;
  int ny = gdcurv->ny;
  int nz = gdcurv->nz;
  int gni1 = gdcurv->gni1;
  int gnj1 = gdcurv->gnj1;
  int gnk1 = gdcurv->gnk1;
  int total_nx = gdcurv->total_nx;
  int total_ny = gdcurv->total_ny;
  int total_nz = gdcurv->total_nz;
  int gni, gnj, gnk;
  size_t siz_iy = gdcurv->siz_iy;
  size_t siz_iz = gdcurv->siz_iz;

  float *P_x1 = src->P_x1;
  float *Q_x1 = src->Q_x1;
  float *R_x1 = src->R_x1;
  float *P_x2 = src->P_x2;
  float *Q_x2 = src->Q_x2;
  float *R_x2 = src->R_x2;

  float *P_y1 = src->P_y1;
  float *Q_y1 = src->Q_y1;
  float *R_y1 = src->R_y1;
  float *P_y2 = src->P_y2;
  float *Q_y2 = src->Q_y2;
  float *R_y2 = src->R_y2;
  
  float *P_z1 = src->P_z1;
  float *Q_z1 = src->Q_z1;
  float *R_z1 = src->R_z1;
  float *P_z2 = src->P_z2;
  float *Q_z2 = src->Q_z2;
  float *R_z2 = src->R_z2;

  float *P = src->P;
  float *Q = src->Q;
  float *R = src->R;

  float xi,et,zt,c0,c1,r0,r1;
  size_t iptr, iptr1;

  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        gni = gni1 + i;
        gnk = gnk1 + k-1; 
        gnj = gnj1 + j-1; 
        xi = (1.0*gni)/(total_nx-1);

        c0 = 1-xi;
        c1 = xi;

        r0 = exp(-coef[0]*xi);
        r1 = exp(-coef[1]*(1-xi)); 
        
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = gnk*(total_ny-2) + gnj;
        P[iptr] = c0*P_x1[iptr1] + c1*P_x2[iptr1];
        Q[iptr] = r0*Q_x1[iptr1] + r1*Q_x2[iptr1];
        R[iptr] = r0*R_x1[iptr1] + r1*R_x2[iptr1];
      }
    }
  }

  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        gnj = gnj1 + j;
        gnk = gnk1 + k-1; 
        gni = gni1 + i-1; 
        et = (1.0*gnj)/(total_ny-1);

        c0 = 1-et;
        c1 = et;

        r0 = exp(-coef[2]*et);
        r1 = exp(-coef[3]*(1-et)); 
        
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = gnk*(total_nx-2) + gni;
        P[iptr] = P[iptr] + r0*P_y1[iptr1] + r1*P_y2[iptr1];
        Q[iptr] = Q[iptr] + c0*Q_y1[iptr1] + c1*Q_y2[iptr1];
        R[iptr] = R[iptr] + r0*R_y1[iptr1] + r1*R_y2[iptr1];
      }
    }
  }


  for(int k=1; k<nz-1; k++) {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i++)
      {
        gnk = gnk1 + k; // include bdry point 
        gnj = gnj1 + j-1; // not include bdry point
        gni = gni1 + i-1;
        zt = (1.0*gnk)/(total_nz-1);
        c0 = 1-zt;
        c1 = zt;

        r0 = exp(-coef[4]*zt);
        r1 = exp(-coef[5]*(1-zt)); 
        
        iptr  = k*siz_iz + j*siz_iy + i;
        iptr1 = gnj*(total_nx-2) + gni;
        P[iptr] = P[iptr] + r0*P_z1[iptr1] + r1*P_z2[iptr1];
        Q[iptr] = Q[iptr] + r0*Q_z1[iptr1] + r1*Q_z2[iptr1];
        R[iptr] = R[iptr] + c0*R_z1[iptr1] + c1*R_z2[iptr1];
      }
    }
  }

  return 0;
}

