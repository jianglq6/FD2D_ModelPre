
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* LOC */
int media_parameterization_loc(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu, float *lam, float *muxz, float *rho_x, float *rho_z);

int media_parameterization_gri(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu, float *lam, float *muxz,  float *Bx, float *Bz, float *rho, float *mu);

int media_parameterization_vol(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu11, float *lam11, float *mu00, float *rho01, float *rho10,
        float *lam2mu, float *lam, float *muxz,  float *rho_x, float *rho_z);

int media_parameterization_tti(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu00, float *lam2mu11, float *lam00, float *lam11, float *mu00, float *mu11, float *rho01, float *rho10,
        float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
        float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
        float *rho_x, float *rho_z);

/* Calculate volume arithmetic averaging */
float cal_ari_averaging(int loc_type, float val00, float val01, float val10, float val11, float area1, float area2);

/* Calculate volume harmonic averaging */
float cal_har_averaging(int loc_type, float val00, float val01, float val10, float val11, float area1, float area2);
