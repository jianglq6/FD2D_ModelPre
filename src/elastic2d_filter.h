#include <stdio.h>
#include <stdlib.h>


/* Selective Filter */
#define SIGMAD 0.2

int filter_velocity(int half_fd_stencil, int nx, int nz,
                    float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2);

int filter_stresses(int half_fd_stencil, int nx, int nz,
                    float *Txx_1, float *Tzz_1, float *Txz_1,
                    float *Txx_2, float *Tzz_2, float *Txz_2);

/* Calculate the filter coefficient, easy to change other filter */
int filter_coe(int half_fd_stencil, float *f_Coe);

/* 2D selective filtering for (i, j) points */
int sf00(int half_sf_stencil, float *sf, float *u00, float *u11, int nx, int nz);

/* 2D selective filtering for (i+1/2, j+1/2) points */
int sf11(int half_sf_stencil, float *sf, float *u00, float *u11, int nx, int nz);

/* 2D selective filtering for (i+1/2, j) points */
int sf01(int half_sf_stencil, float *sf, float *u01, float *u10, int nx, int nz);

/* 2D selective filtering for (i, j+1/2) points */
int sf10(int half_sf_stencil, float *sf, float *u01, float *u10, int nx, int nz);

int SFsCoe(int half_fd_stencil, float *SFs);

int SFoCoe(int half_fd_stencil, float *SFo);
