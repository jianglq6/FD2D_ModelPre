#include <stdio.h>
#include <stdlib.h>


#define NOFILTER   0
#define FILTER_SFs 1
#define FILTER_SFo 2

/* Selective Filter */
#define SIGMAD 0.2

int filter_velocity(int filter_method, int half_fd_stencil, int nx,
                    int ni1, int ni2, int nk1, int nk2,
                    float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2);

int filter_stresses(int filter_method, int half_fd_stencil, int nx,
                    int ni1, int ni2, int nk1, int nk2,
                    float *Txx_1, float *Tzz_1, float *Txz_1,
                    float *Txx_2, float *Tzz_2, float *Txz_2);

/* Calculate the filter coefficient, easy to change other filter */
int filter_coe(int filter_method, int half_fd_stencil, float *f_Coe, int *half_sf_stencil);

/* 2D selective filtering for (i, j) points */
int sf00(int half_sf_stencil, float *sf, float *u00, float *u11,
         int nx, int ni1, int ni2, int nk1, int nk2);

/* 2D selective filtering for (i+1/2, j+1/2) points */
int sf11(int half_sf_stencil, float *sf, float *u00, float *u11,
         int nx, int ni1, int ni2, int nk1, int nk2);

/* 2D selective filtering for (i+1/2, j) points */
int sf01(int half_sf_stencil, float *sf, float *u01, float *u10,
         int nx, int ni1, int ni2, int nk1, int nk2);

/* 2D selective filtering for (i, j+1/2) points */
int sf10(int half_sf_stencil, float *sf, float *u01, float *u10,
         int nx, int ni1, int ni2, int nk1, int nk2);

int SFsCoe(int half_fd_stencil, float *SFs, int *half_sf_stencil);

int SFoCoe(int half_fd_stencil, float *SFo, int *half_sf_stencil);
