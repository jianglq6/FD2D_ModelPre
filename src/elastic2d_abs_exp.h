#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ALPHA 0.025
#define MIN(x,y) ((x) < (y) ? (x) : (y))

int abs_exp_velocity(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
                     int *boundary_type, int *boundary_layer_number,
                     float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2);

int abs_exp_stresses(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
                     int *boundary_type, int *boundary_layer_number,
                     float *Txx_1, float *Txx_2, float *Tzz_1, float *Tzz_2, float *Txz_1, float *Txz_2);

int abs_exp_x1(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U);

int abs_exp_x2(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U);

int abs_exp_z1(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U);

int abs_exp_z2(int ni1, int ni2, int nk1, int nk2, int nx, int nz, int half_fd_stencil,
               int boundary_layer_number, float *U);
