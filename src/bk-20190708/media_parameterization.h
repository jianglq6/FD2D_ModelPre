
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* LOC */
int media_parameterization_loc(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *xvec, float *zvec,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam, float *mu, float *muxz, float *rho_x, float *rho_z);

int media_parameterization_gri(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *xvec, float *zvec,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam, float *lam2mu, float *muxz,  float *Bx, float *Bz, float *rho, float *mu);

int media_parameterization_vol(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *xvec, float *zvec,
        float *LAM, float *MU, float *RHO, float *lam11, float *mu00, float *mu11, float *rho01, float *rho10,
        float *lam, float *mu, float *muxz,  float *rho_x, float *rho_z);

float cal_para_vol(int loc_type, float val00, float val01, float val10, float val11,
    float area1, float area2);
