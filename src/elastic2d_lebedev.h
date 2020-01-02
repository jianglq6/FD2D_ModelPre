#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


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

//float Dx00( float *a, double *fdx, int ith, int indx, int nx);
//float Dx01( float *a, double *fdx, int ith, int indx, int nx);
//float Dx10( float *a, double *fdx, int ith, int indx, int nx);
//float Dx11( float *a, double *fdx, int ith, int indx, int nx);
//
//float Dz00( float *a, double *fdz, int ith, int indx, int nx);
//float Dz01( float *a, double *fdz, int ith, int indx, int nx);
//float Dz10( float *a, double *fdz, int ith, int indx, int nx);
//float Dz11( float *a, double *fdz, int ith, int indx, int nx);


int elastic2d_lebedev(float dx, float dz, int nx, int nz, int nt, float dt,
                    int half_fd_stencil, int spatial_difference_method, int filter_method,
                    float *xvec1, float *zvec1, float *xvec2, float *zvec2,
                    float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
                    float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
                    float *B01, float *B10, struct Src src, int source_impulse_method,
                    int abs_type, int *boundary_layer_number,
                    float *Txx00, float *Txx11, float *Txz11, float *Txz00, float *Tzz00, float *Tzz11,
                    float *Vx01, float *Vx10, float *Vz10, float *Vz01,
                    float *hTxx00, float *hTxx11, float *hTxz11, float *hTxz00, float *hTzz00, float *hTzz11,
                    float *hVx01, float *hVx10, float *hVz10, float *hVz01,
                    float *DxVx00, float *DzVz00, float *DzVx00, float *DxVz00,
                    float *DxVx11, float *DzVz11, float *DzVx11, float *DxVz11,
                    float *DxTxx01,float *DzTxz01,float *DxTxz01,float *DzTzz01,
                    float *DxTxx10,float *DzTxz10,float *DxTxz10,float *DzTzz10,
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
int cal_momentum(double *fdx, double *fdz, int half_fd_stencil,
                 int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx01, float *hVz10, float *Txx00, float *Txz11, float *Tzz00,
                 float *hVx10, float *hVz01, float *Txx11, float *Txz00, float *Tzz11,
                 float *DxTxx01, float *DzTxz01, float *DxTxz01, float *DzTzz01,
                 float *DxTxx10, float *DzTxz10, float *DxTxz10, float *DzTzz10);

/* RHS of stress-strain equation (Hook law) */
int cal_hook(double *fdx, double *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx00, float *hTzz00, float *hTxz11, float *Vx01 , float *Vz10,
             float *hTxx11, float *hTzz11, float *hTxz00, float *Vx10 , float *Vz01,
             float *DxVx00, float *DzVz00, float *DzVx00, float *DxVz00,
             float *DxVx11, float *DzVz11, float *DzVx11, float *DxVz11,
             float *c11_1 , float *c13_1 , float *c15_1 , float *c33_1, float *c35_1, float *c55_1,
             float *c11_2 , float *c13_2 , float *c15_2 , float *c33_2, float *c35_2, float *c55_2);

int update_velocity(int nx, int nz, float dt,
        float *Vx01 , float *Vz10 , float *Vx10 , float *Vz01 ,
        float *hVx01, float *hVz10, float *hVx10, float *hVz01,
        float *B10, float *B01);

int update_stress(int nx, int nz, float dt,
                  float *Txx00 , float *Tzz00 , float *Txz11,
                  float *Txx11 , float *Tzz11 , float *Txz00,
                  float *hTxx00, float *hTzz00, float *hTxz11,
                  float *hTxx11, float *hTzz11, float *hTxz00 );
