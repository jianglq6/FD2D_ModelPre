#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define BOUNDARY_TYPE_FREESURFACE 0
#define BOUNDARY_TYPE_EXPONENT 1

/* ith-difference to 00 position  */
#define Dx00(a,fdx,i_diff,ix,iz,nx)   \
    ( fdx[i_diff-1] * ( a[iz*nx + ix+(i_diff-1)] - a[iz*nx + ix-(i_diff  ) ] ) )

#define Dx10(a,fdx,i_diff,ix,iz,nx)   \
    ( fdx[i_diff-1] * ( a[iz*nx + ix+(i_diff-1)] - a[iz*nx + ix-(i_diff  ) ] ) )

#define Dx11(a,fdx,i_diff,ix,iz,nx)   \
    ( fdx[i_diff-1] * ( a[iz*nx + ix+(i_diff  )] - a[iz*nx + ix-(i_diff-1) ] ) )

#define Dx01(a,fdx,i_diff,ix,iz,nx)   \
    ( fdx[i_diff-1] * ( a[iz*nx + ix+(i_diff  )] - a[iz*nx + ix-(i_diff-1) ] ) )

#define Dz00(a,fdz,i_diff,ix,iz,nx)   \
    ( fdz[i_diff-1] * ( a[(iz+(i_diff-1))*nx + ix] - a[(iz- i_diff   )*nx + ix] ) )

#define Dz01(a,fdz,i_diff,ix,iz,nx)   \
    ( fdz[i_diff-1] * ( a[(iz+(i_diff-1))*nx + ix] - a[(iz- i_diff   )*nx + ix] ) )

#define Dz11(a,fdz,i_diff,ix,iz,nx)   \
    ( fdz[i_diff-1] * ( a[(iz+ i_diff   )*nx + ix] - a[(iz-(i_diff-1))*nx + ix] ) )

#define Dz10(a,fdz,i_diff,ix,iz,nx)   \
    ( fdz[i_diff-1] * ( a[(iz+ i_diff   )*nx + ix] - a[(iz-(i_diff-1))*nx + ix] ) )

int elastic2d_staggered(float dx, float dz, int nx, int nz, int nt, float dt,
                    int half_fd_stencil, int spatial_difference_method,
                    float *xvec1, float *zvec1, float *xvec2, float *zvec2,
                    float *c11, float *c13, float *c33, float *c55,
                    float *B01, float *B10,
                    struct Src src, int source_impulse_method,
                    int abs_type, int *boundary_layer_number,
                    float *Txx, float *Txz, float *Tzz,
                    float *Vx, float *Vz,
                    float *hTxx, float *hTxz, float *hTzz,
                    float *hVx, float *hVz,
                    float *DxVx00, float *DzVz00, float *DzVx11, float *DxVz11,
                    float *DxTxx01,float *DzTxz01,float *DxTxz10,float *DzTzz10,
                    int seismotype, int nreceivers, float *xr, float *zr,
                    int NSTEP_BETWEEN_OUTPUT_SEISMOS,
                    bool save_ASCII_seismograms, bool save_binary_seismograms,
                    int NSTEP_BETWEEN_OUTPUT_INFO,
                    int NSTEP_BETWEEN_OUTPUT_IMAGES,
                    bool output_postscript_snapshot,
                    int imagetype_postscript, bool meshvect, bool modelvect, bool boundvect,
                    float sizemax_arrows, bool US_LETTER,
                    bool output_wavefield_dumps, int imagetype_wavefield_dumps,
                    bool use_binary_for_wavefield_dumps);

/* RHS of equation of motion */
int cal_momentum_ssg(double *fdx, double *fdz, int half_fd_stencil,
                 int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx, float *hVz, float *Txx, float *Txz, float *Tzz,
                 float *DxTxx01, float *DzTxz01, float *DxTxz10, float *DzTzz10);

/* RHS of stress-strain equation (Hook law) */
int cal_hook_ssg(double *fdx, double *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx, float *hTzz, float *hTxz, float *Vx , float *Vz,
             float *c11 , float *c13 , float *c33 , float *c55,
             float *DxVx00, float *DzVz00, float *DxVz11, float *DzVx11);

int update_velocity_ssg(int nx, int nz, float dt,
        float *Vx, float *Vz,
        float *hVx, float *hVz,
        float *B10, float *B01);

int update_stress_ssg(int nx, int nz, float dt,
                  float *Txx , float *Tzz , float *Txz,
                  float *hTxx, float *hTzz, float *hTxz);

