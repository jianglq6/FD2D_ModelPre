#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ABS_NPML 2

int abs_npml_velocity_tti(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                          int *abs_layer_number, float dt,
                          float *hVx01,float *hVx10,float *hVz10,float *hVz01,
                          float *Ax_regu,  float *Bx_regu, float *Dx_regu,
                          float *Ax_half,  float *Bx_half, float *Dx_half,
                          float *Az_regu,  float *Bz_regu, float *Dz_regu,
                          float *Az_half,  float *Bz_half, float *Dz_half,
                          float *Txx01_x1, float *Txx01_x2,
                          float *Tzz01_z1, float *Tzz01_z2,
                          float *Txz01_x1, float *Txz01_x2,
                          float *Txz01_z1, float *Txz01_z2,
                          float *Txx10_x1, float *Txx10_x2,
                          float *Tzz10_z1, float *Tzz10_z2,
                          float *Txz10_x1, float *Txz10_x2,
                          float *Txz10_z1, float *Txz10_z2,
                          float *DxTxx01,  float *DzTxz01, float *DxTxz01, float *DzTzz01,
                          float *DxTxx10,  float *DzTxz10, float *DxTxz10, float *DzTzz10);

int abs_npml_stress_tti(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                        int *abs_layer_number, float dt,
                        float *hTxx00,  float *hTxx11,
                        float *hTxz00,  float *hTxz11,
                        float *hTzz00,  float *hTzz11,
                        float *DxVx00,  float *DzVz00,  float *DzVx00,  float *DxVz00,
                        float *DxVx11,  float *DzVz11,  float *DzVx11,  float *DxVz11,
                        float *c11_1 ,  float *c13_1 ,  float *c15_1 ,  float *c33_1, float *c35_1, float *c55_1,
                        float *c11_2 ,  float *c13_2 ,  float *c15_2 ,  float *c33_2, float *c35_2, float *c55_2,
                        float *Ax_regu, float *Bx_regu, float *Dx_regu,
                        float *Ax_half, float *Bx_half, float *Dx_half,
                        float *Az_regu, float *Bz_regu, float *Dz_regu,
                        float *Az_half, float *Bz_half, float *Dz_half,
                        float *Vx00_x1, float *Vx00_x2, float *Vx11_x1, float *Vx11_x2,
                        float *Vx00_z1, float *Vx00_z2, float *Vx11_z1, float *Vx11_z2,
                        float *Vz00_x1, float *Vz00_x2, float *Vz11_x1, float *Vz11_x2,
                        float *Vz00_z1, float *Vz00_z2, float *Vz11_z1, float *Vz11_z2);

int abs_npml_init(int *abs_layer_number,
                  float fc, float *c11, float *c55, float *rho_recip,
                  int half_fd_stencil,
                  float *Ax_regu, float *Bx_regu, float *Dx_regu,
                  float *Ax_half, float *Bx_half, float *Dx_half,
                  float *Az_regu, float *Bz_regu, float *Dz_regu,
                  float *Az_half, float *Bz_half, float *Dz_half,
                  float *xvec1, float *xvec2, float *zvec1, float *zvec2,
                  int ni1, int ni2, int nk1, int nk2, int nx, float dx, float dz);

float cal_pml_R(int N);

float cal_pml_d0(float L, float vp, float R);

float cal_pml_alpha0(float fc);

float cal_pml_beta0(float PPW_s, float PPW_FD);

float cal_pml_d(float d0, float x, float L);

float cal_pml_alpha(float alpha0, float x, float L);

float cal_pml_beta(float beta0, float x, float L);

int cal_vel_pml(float *c11, float *rho_recip, float *vp_pml,
               int *abs_layer_number,
               int ni1, int ni2, int nk1, int nk2, int nx);

int abs_npml_velocity_staggered(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                                int *abs_layer_number, float dt,
                                float *hVx, float *hVz,
                                float *Ax_regu, float *Bx_regu, float *Dx_regu,
                                float *Ax_half, float *Bx_half, float *Dx_half,
                                float *Az_regu, float *Bz_regu, float *Dz_regu,
                                float *Az_half, float *Bz_half, float *Dz_half,
                                float *Txx_x1, float *Txx_x2,
                                float *Tzz_z1, float *Tzz_z2,
                                float *Txz_x1, float *Txz_x2,
                                float *Txz_z1, float *Txz_z2,
                                float *DxTxx01, float *DzTxz01,
                                float *DxTxz10, float *DzTzz10);

int abs_npml_stress_staggered(int ni1, int ni2, int nk1, int nk2, int nx, int nz,
                              int *abs_layer_number, float dt,
                              float *hTxx, float *hTxz, float *hTzz,
                              float *DxVx00, float *DzVz00, float *DzVx11, float *DxVz11,
                              float *c11, float *c13, float *c33, float *c55,
                              float *Ax_regu, float *Bx_regu, float *Dx_regu,
                              float *Ax_half, float *Bx_half, float *Dx_half,
                              float *Az_regu, float *Bz_regu, float *Dz_regu,
                              float *Az_half, float *Bz_half, float *Dz_half,
                              float *Vx_x1, float *Vx_x2, float *Vx_z1, float *Vx_z2,
                              float *Vz_x1, float *Vz_x2, float *Vz_z1, float *Vz_z2);


