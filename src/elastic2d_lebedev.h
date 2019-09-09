#include <stdio.h>
#include <stdlib.h>

#define BOUNDARY_TYPE_FREESURFACE 0
#define BOUNDARY_TYPE_EXPONENT 1

/* ith-difference to 00 position  */
#define Dx00(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith-1)]   - a[iz*nx + ix-(ith  ) ] ) )

#define Dx10(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith-1)]   - a[iz*nx + ix-(ith  ) ] ) )

#define Dx11(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith  )]   - a[iz*nx + ix-(ith-1) ] ) )

#define Dx01(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith  )]   - a[iz*nx + ix-(ith-1) ] ) )

#define Dz00(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+(ith-1))*nx + ix] - a[(iz- ith   )*nx + ix] ) )

#define Dz01(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+(ith-1))*nx + ix] - a[(iz- ith   )*nx + ix] ) )

#define Dz11(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+ ith   )*nx + ix] - a[(iz-(ith-1))*nx + ix] ) )

#define Dz10(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+ ith   )*nx + ix] - a[(iz-(ith-1))*nx + ix] ) )


int cal_momentum(float *fdx, float *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx_1, float *hVz_1, float *Txx_1, float *Txz_1, float *Tzz_1,
                 float *hVx_2, float *hVz_2, float *Txx_2, float *Txz_2, float *Tzz_2);

/* RHS of stress-strain equation (Hook law) */
int cal_hook(float *fdx, float *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx_1, float *hTzz_1, float *hTxz_1, float *Vx_1 , float *Vz_1,
             float *hTxx_2, float *hTzz_2, float *hTxz_2, float *Vx_2 , float *Vz_2,
             float *c11_1 , float *c13_1 , float *c15_1 , float *c33_1, float *c35_1, float *c55_1,
             float *c11_2 , float *c13_2 , float *c15_2 , float *c33_2, float *c35_2, float *c55_2);

int update_velocity(int nx, int nz, float dt, float *Vx_1, float *Vz_1, float *Vx_2, float *Vz_2,
        float *hVx_1, float *hVz_1,  float *hVx_2, float *hVz_2, float *B10, float *B01);

int update_stress(int nx, int nz, float dt, float *Txx_1, float *Tzz_1, float *Txz_1,
                  float *Txx_2 , float *Tzz_2 , float *Txz_2 ,
                  float *hTxx_1, float *hTzz_1, float *hTxz_1,
                  float *hTxx_2, float *hTzz_2, float *hTxz_2 );
