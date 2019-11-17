/*************************************************************************
 *
 * These functions is used for media parameterization, which include many
 *  difference method:
 * 1. GRI: grid average method (Graves - 1996 - Simulating seismic wave
 *         propagation in 3D elastic media using staggered-grid finite
 *         differences)
 * 2. VOL: volume arithmetic and harmonic averaging method (Moczo et al.,
 *         - 2002 - 3D heterogeneous staggered-grid finite-difference
 *         modeling of seismic motion with volume harmonic and arithmetic
 *         averaging of elastic moduli and densities)
 * 3. ORT: effective orthorhombic averaged medium method (Kristek et al.,
 *         - 2017 - An orthorhombic representation of a heterogeneous medium
 *         for the finite-difference modelling of seismic wave propagation)
 * 4. TTI: TTI effective media parameterization method
 *
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018-2019 zwlab
 *
 * Version: 1.0
 *
 * Date: 04/2019
 *
 * History:
 *     04/2019: Original version created by Luqian Jiang
 *     11/2019: Luqian Jiang:
 *         combine with modeling, bounduary expansion
 *
 ************************************************************************/

#include "pre_interface_struct.h"
#include "pre_cal_dip_area.h"
#include "pre_media_parameterization.h"
#include "pre_para_assign.h"

/*********************************************************************
 * For LOC, VOL, GRI:
 *   lam, lam2mu - (j, i)
 *   rho_x       - (j, i+1/2)
 *   rho_z       - (j+1/2, i)
 *   mu_xz       - (j+1/2, i+1/2)
 ********************************************************************/

/* LOC */
int media_parameterization_loc(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu, float *lam, float *muxz, float *rho_x, float *rho_z,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    printf("LOC: dx: %f, dz: %f, num_of_interfaces: %d\n",dx,dz,number_of_interfaces);

    int ierr=0, ix, iz;

    // Need error judgement !!!
    /*lam2mu, lam*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0],
            dx, dz, nx, nz, LAM, lam, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0],
            dx, dz, nx, nz, L2M, lam2mu, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /*rho_x*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec1[0],
            dx, dz, nx, nz, RHO, rho_x, nghost_x1, nghost_x2, nghost_z1, nghost_z2);


    /*rho_z*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec2[0],
            dx, dz, nx, nz, RHO, rho_z, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /*muxz*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0],
            dx, dz, nx, nz, MU, muxz, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    return ierr;
}


/* GRI*/
int media_parameterization_gri(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu, float *lam, float *muxz,  float *Bx, float *Bz, float *rho, float *mu,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    printf("GRI: dx: %f, dz: %f, num_of_interfaces: %d\n",dx,dz,number_of_interfaces);

    int ix, iz;
    int ierr=0;


    /*lam2mu, lam*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, L2M, lam2mu, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, LAM, lam, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /* rho */
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, RHO, rho, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    if (!ierr) {
        /* Bx */
        for (ix = nghost_x1; ix < nx-1-nghost_x2; ix++) {
            for (iz = nghost_z1; iz < nz-nghost_z2; iz++) {
                Bx[iz*nx + ix] = 1/rho[iz*nx + ix] + 1/rho[iz*nx + (ix+1)];
            }
        }
        for (iz = nghost_z1; iz < nz-nghost_z2; iz++) {
            Bx[iz*nx + nx] = Bx[iz*nx + (nx-1)];
        }

        /* Bz */
        for (ix = nghost_x1; ix < nx-nghost_x2; ix++) {
            for (iz = nghost_z1+1; iz < nz-nghost_z2; iz++) {
                Bz[iz*nx + ix] = 1/rho[iz*nx + ix] + 1/rho[(iz-1)*nx + ix];
            }
        }
        for (ix = nghost_x1; ix < nx-nghost_x2; ix++) {
            Bz[0*nx + ix] = Bz[1*nx + ix];
        }
    }

    /* muxz */
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, MU, mu, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    for (ix = nghost_x1; ix < nx-1-nghost_x2; ix++) {
        for (iz = nghost_z1+1; iz < nz-nghost_z2; iz++) {
            muxz[iz*nx + ix] = 4.0 / ( 1/mu[(iz-1)*nx + ix] + 1/mu[iz*nx + ix    ]
                                    +  1/mu[(iz-1)*nx + ix] + 1/mu[iz*nx + (ix+1)]);
        }
    }
    for (iz = nghost_z1 + 1; iz < nz-nghost_z2; iz++) {
        muxz[iz*nx + nx] = mu[iz*nx + (nx-1)];
    }
    for (ix = nghost_x1; ix < nx - nghost_x2; ix++) {
        muxz[0*nx + ix] = mu[1*nx + ix];
    }


    free(rho);
    free(mu);

    return ierr;

}


int media_parameterization_vol(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu11, float *lam11, float *mu00, float *rho01, float *rho10,
        float *lam2mu, float *lam, float *muxz,  float *rho_x, float *rho_z,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    int   *ix = NULL, *iz = NULL, *loc_type = NULL; //ix, iz there is the layer of grid
    float *area1 = NULL, *area2 = NULL, *theta = NULL;
    float xvec_h[nx], zvec_h[nz];  // xvec-dx/2, zvec-dz/2.
    int   ix_h, iz_h, npoints_layer, i_point, index, ierr, i_interface;

    int i,j;

    for (ix_h = 0; ix_h < nx; ix_h++)
        xvec_h[ix_h] = xvec1[ix_h] - dx/2;
    for (iz_h = 0; iz_h < nz; iz_h++)
        zvec_h[iz_h] = zvec1[iz_h] - dz/2;

    ix       = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    iz       = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    loc_type = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    area1    = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    area2    = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    theta    = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));

    memset(ix, 0, MAX_LAYER_POINTS*sizeof(int));
    memset(iz, 0, MAX_LAYER_POINTS*sizeof(int));
    memset(loc_type, 0, MAX_LAYER_POINTS*sizeof(int));
    memset(area1, 0.0, MAX_LAYER_POINTS*sizeof(float));
    memset(area2, 0.0, MAX_LAYER_POINTS*sizeof(float));
    memset(theta, 0.0, MAX_LAYER_POINTS*sizeof(float));


    /* assginment first */
    /* mu, lam*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, L2M, lam2mu, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, LAM, lam, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, MU, mu00, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /*rho_x*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec1[0], dx, dz,
            nx, nz, RHO, rho_x, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec1[0], dx, dz,
            nx, nz, RHO, rho01, nghost_x1, nghost_x2, nghost_z1, nghost_z2);


    /*rho_z*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec2[0], dx, dz,
            nx, nz, RHO, rho_z, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec2[0], dx, dz,
            nx, nz, RHO, rho10, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /*muxz*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, MU, muxz, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, L2M, lam2mu11, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, LAM, lam11, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /* volumn arithmetic and harmonic averaging */
    /* calculate the area */
    for (i_interface = 0; i_interface < number_of_interfaces; i_interface++) {
        printf("Interface number: %d\n",i_interface);
        /*-- For sigma_xx, sigma_zz: lam, mu, center of the volum is 00 --*/
        cal_dip_area(xvec_h, zvec_h, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] >= nghost_z1 && iz[i_point] < nz-nghost_z2
              && ix[i_point] >= nghost_x1 && ix[i_point] < nx-nghost_x2 && loc_type[i_point] != 0 ) {
                index = iz[i_point]*nx+ix[i_point];
                lam[index] = cal_har_averaging(loc_type[i_point], lam11[index-1-nx], lam11[index-nx], lam11[index-1],
                                        lam11[index], area1[i_point], area2[i_point]);
                lam2mu[index] = cal_har_averaging(loc_type[i_point], lam2mu11[index-1-nx], lam2mu11[index-nx],
                                        lam2mu11[index-1], lam2mu11[index], area1[i_point], area2[i_point]);

            }

        }

        /*-- For sigma_xz: muxz , center of the volum is 11 --*/
        cal_dip_area(xvec1, zvec1, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] >= nghost_z1 && iz[i_point] < nz-nghost_z2
              && ix[i_point] >= nghost_x1 && ix[i_point] < nx-nghost_x2 && loc_type[i_point] != 0 ) {
                index = iz[i_point]*nx+ix[i_point];
                muxz[index] = cal_har_averaging(loc_type[i_point], mu00[index], mu00[index+1], mu00[index+nx],
                                                mu00[index+nx+1], area1[i_point], area2[i_point]);

            }
        }

      /*-- For v_x: rho_x , center of the volum is 01 --*/
        cal_dip_area(xvec1, zvec_h, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] >= nghost_z1 && iz[i_point] < nz-nghost_z2
              && ix[i_point] >= nghost_x1 && ix[i_point] < nx-nghost_x2 && loc_type[i_point] != 0 ) {
                index = iz[i_point]*nx+ix[i_point];
                rho_x[index] = cal_ari_averaging(loc_type[i_point], rho10[index-nx], rho10[index-nx+1], rho10[index],
                                                 rho10[index+1], area1[i_point], area2[i_point]);
            };
        }

        /*-- For v_z: rho_z , center of the volum is 10 --*/
        cal_dip_area(xvec_h, zvec1, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer ,ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] >= nghost_z1 && iz[i_point] < nz-nghost_z2
              && ix[i_point] >= nghost_x1 && ix[i_point] < nx-nghost_x2 && loc_type[i_point] != 0 ) {
                index = iz[i_point]*nx+ix[i_point];
                rho_z[index] = cal_ari_averaging(loc_type[i_point], rho01[index-1], rho01[index], rho01[index+nx-1],
                                                 rho01[index+nx], area1[i_point], area2[i_point]);
            }
        }

    }

    free(ix);
    free(iz);
    free(loc_type);
    free(area1);
    free(area2);
    free(theta);
    return 0;
}

int media_parameterization_tti(struct interfaces *interface, int number_of_interfaces,
        float dx, float dz, int nx, int nz, float *xvec1, float *zvec1, float *xvec2, float *zvec2,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam2mu00, float *lam2mu11, float *lam00, float *lam11, float *mu00, float *mu11, float *rho01, float *rho10,
        float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
        float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
        float *rho_x, float *rho_z,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    int   *ix = NULL, *iz = NULL, *loc_type = NULL;  //ix, iz there is the layer of grid
    float *area1 = NULL, *area2 = NULL, *theta = NULL;
    float xvec_h[nx], zvec_h[nz];  // xvec-dx/2, zvec-dz/2.
    int   ix_h, iz_h, npoints_layer, i_point, index, ierr, i_interface;
    float A1, B1, C1, A2, B2, C2, tmp1, tmp2, muxz00, muxz11;
    int i,j;

    for (ix_h = 0; ix_h < nx; ix_h++)
        xvec_h[ix_h] = xvec1[ix_h] - dx/2;
    for (iz_h = 0; iz_h < nz; iz_h++)
        zvec_h[iz_h] = zvec1[iz_h] - dz/2;

    ix       = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    iz       = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    loc_type = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    area1    = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    area2    = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    theta    = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));

    memset(ix, 0, MAX_LAYER_POINTS*sizeof(int));
    memset(iz, 0, MAX_LAYER_POINTS*sizeof(int));
    memset(loc_type, 0, MAX_LAYER_POINTS*sizeof(int));
    memset(area1, 0, MAX_LAYER_POINTS*sizeof(float));
    memset(area2, 0, MAX_LAYER_POINTS*sizeof(float));
    memset(theta, 0.0, MAX_LAYER_POINTS*sizeof(float));


    /* assginment all the grid */
    /* mu, lam*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, MU , mu00, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, LAM, lam00, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, L2M, lam2mu00, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, L2M, c11_1, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, L2M, c33_1, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, LAM, c13_1, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec1[0], dx, dz,
            nx, nz, MU , c55_2, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    memset(c15_1, 0.0, nx*nz*sizeof(float));
    memset(c35_1, 0.0, nx*nz*sizeof(float));

    /*rho_x*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec1[0], dx, dz,
            nx, nz, RHO, rho_x, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec1[0], dx, dz,
            nx, nz, RHO, rho01, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /*rho_z*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec2[0], dx, dz,
            nx, nz, RHO, rho_z, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec1[0], zvec2[0], dx, dz,
            nx, nz, RHO, rho10, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    /*muxz*/
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, MU , mu11, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, LAM, lam11, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, L2M, lam2mu11, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, L2M, c11_2, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, L2M, c33_2, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, LAM, c13_2, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    ierr = parameter_assignment(interface, number_of_interfaces, xvec2[0], zvec2[0], dx, dz,
            nx, nz, MU , c55_1, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    memset(c15_2, 0.0, nx*nz*sizeof(float));
    memset(c35_2, 0.0, nx*nz*sizeof(float));


    /* TTI Approximation */
    /* calculate the area */
    for (i_interface = 0; i_interface < number_of_interfaces; i_interface++) {
        printf("Interface number: %d\n",i_interface);
        /*-- For sigma_xx, sigma_zz: lam, mu, center of the volum is 00 --*/
        cal_dip_area(xvec_h, zvec_h, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer,ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > nghost_z1 && iz[i_point] < nz-nghost_z2-1
              && ix[i_point] > nghost_x1 && ix[i_point] < nx-nghost_x2-1 && loc_type[i_point] != 0 ) {

                index = iz[i_point]*nx+ix[i_point];

                muxz00 = cal_har_averaging(loc_type[i_point], mu11[index-1-nx], mu11[index-nx], mu11[index-1],
                                           mu11[index], area1[i_point], area2[i_point]);

                tmp1   = cal_ari_averaging(loc_type[i_point], lam11[index-1-nx]/lam2mu11[index-1-nx],
                                           lam11[index-nx]/lam2mu11[index-nx], lam11[index-1]/lam2mu11[index-1],
                                           lam11[index]/lam2mu11[index], area1[i_point], area2[i_point]);

                tmp2   = cal_ari_averaging(loc_type[i_point],
                                           lam2mu11[index-1-nx] - pow(lam11[index-1-nx], 2)/lam2mu11[index-1-nx],
                                           lam2mu11[index-nx  ] - pow(lam11[index-nx  ], 2)/lam2mu11[index-nx  ],
                                           lam2mu11[index-1   ] - pow(lam11[index-1   ], 2)/lam2mu11[index-1   ],
                                           lam2mu11[index     ] - pow(lam11[index     ], 2)/lam2mu11[index     ],
                                           area1[i_point], area2[i_point]);

                A1     = cal_har_averaging(loc_type[i_point], lam2mu11[index-1-nx], lam2mu11[index-nx], lam2mu11[index-1],
                                           lam2mu11[index], area1[i_point], area2[i_point]);

                B1     = tmp1 * A1;
                C1     = tmp2 + tmp1*tmp1 * A1;

                c11_1[index] =  ( C1 * pow(cos(theta[i_point]),2) + B1 * pow(sin(theta[i_point]),2) ) * pow(cos(theta[i_point]),2)
                              + ( B1 * pow(cos(theta[i_point]),2) + A1 * pow(sin(theta[i_point]),2) ) * pow(sin(theta[i_point]),2)
                              + muxz00 * pow(sin(2*theta[i_point]),2);
                c13_1[index] =  ( C1 * pow(cos(theta[i_point]),2) + B1 * pow(sin(theta[i_point]),2) ) * pow(sin(theta[i_point]),2)
                              + ( B1 * pow(cos(theta[i_point]),2) + A1 * pow(sin(theta[i_point]),2) ) * pow(cos(theta[i_point]),2)
                              - muxz00 * pow(sin(2*theta[i_point]),2);
                c33_1[index] =  ( A1 * pow(cos(theta[i_point]),2) + B1 * pow(sin(theta[i_point]),2) ) * pow(cos(theta[i_point]),2)
                              + ( B1 * pow(cos(theta[i_point]),2) + C1 * pow(sin(theta[i_point]),2) ) * pow(sin(theta[i_point]),2)
                              + muxz00 * pow(sin(2*theta[i_point]),2);
                c15_1[index] =  ( C1 * pow(cos(theta[i_point]),2) + B1 * pow(sin(theta[i_point]),2)
                                - B1 * pow(cos(theta[i_point]),2) - A1 * pow(sin(theta[i_point]),2) ) * sin(2*theta[i_point])/2
                              - muxz00 * sin(4*theta[i_point])/2;
                c35_1[index] =  ( B1 * pow(cos(theta[i_point]),2) + C1 * pow(sin(theta[i_point]),2)
                                - A1 * pow(cos(theta[i_point]),2) - B1 * pow(sin(theta[i_point]),2) ) * sin(2*theta[i_point])/2
                              + muxz00 * sin(4*theta[i_point])/2;
                c55_2[index] =  ( C1 - 2*B1 + A1) * pow(sin(2*theta[i_point]),2)/4 + muxz00 * pow(cos(2*theta[i_point]),2);

            }

        }

        /*-- For sigma_xz: muxz , center of the volum is 11 --*/
        cal_dip_area(xvec1, zvec1, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > nghost_z1 && iz[i_point] < nz-nghost_z2-1
              && ix[i_point] > nghost_x1 && ix[i_point] < nx-nghost_x2-1 && loc_type[i_point] != 0 ) {

                index = iz[i_point]*nx+ix[i_point];

                muxz11 = cal_har_averaging(loc_type[i_point], mu00[index], mu00[index+1], mu00[index+nx],
                                                  mu00[index+nx+1], area1[i_point], area2[i_point]);

                tmp1   = cal_ari_averaging(loc_type[i_point], lam00[index]/lam2mu00[index],
                                           lam00[index+1]/lam2mu00[index+1], lam00[index+nx]/lam2mu00[index+nx],
                                           lam00[index+nx+1]/lam2mu00[index+nx+1], area1[i_point], area2[i_point]);

                tmp2   = cal_ari_averaging(loc_type[i_point],
                                           lam2mu00[index     ] - pow(lam00[index     ], 2)/lam2mu00[index     ],
                                           lam2mu00[index+1   ] - pow(lam00[index+1   ], 2)/lam2mu00[index+1   ],
                                           lam2mu00[index+nx  ] - pow(lam00[index+nx  ], 2)/lam2mu00[index+nx  ],
                                           lam2mu00[index+nx+1] - pow(lam00[index+nx+1], 2)/lam2mu00[index+nx+1],
                                           area1[i_point], area2[i_point]);

                A2     = cal_har_averaging(loc_type[i_point], lam2mu00[index], lam2mu00[index+1], lam2mu00[index+nx],
                                           lam2mu00[index+nx+1], area1[i_point], area2[i_point]);

                B2     = tmp1 * A2;
                C2     = tmp2 + tmp1*tmp1 * A2;

                c11_2[index] =  ( C2 * pow(cos(theta[i_point]),2) + B2 * pow(sin(theta[i_point]),2) ) * pow(cos(theta[i_point]),2)
                              + ( B2 * pow(cos(theta[i_point]),2) + A2 * pow(sin(theta[i_point]),2) ) * pow(sin(theta[i_point]),2)
                              + muxz11 * pow(sin(2*theta[i_point]),2);
                c13_2[index] =  ( C2 * pow(cos(theta[i_point]),2) + B2 * pow(sin(theta[i_point]),2) ) * pow(sin(theta[i_point]),2)
                              + ( B2 * pow(cos(theta[i_point]),2) + A2 * pow(sin(theta[i_point]),2) ) * pow(cos(theta[i_point]),2)
                              - muxz11 * pow(sin(2*theta[i_point]),2);
                c33_2[index] =  ( A2 * pow(cos(theta[i_point]),2) + B2 * pow(sin(theta[i_point]),2) ) * pow(cos(theta[i_point]),2)
                              + ( B2 * pow(cos(theta[i_point]),2) + C2 * pow(sin(theta[i_point]),2) ) * pow(sin(theta[i_point]),2)
                              + muxz11 * pow(sin(2*theta[i_point]),2);
                c15_2[index] =  ( C2 * pow(cos(theta[i_point]),2) + B2 * pow(sin(theta[i_point]),2)
                                - B2 * pow(cos(theta[i_point]),2) - A2 * pow(sin(theta[i_point]),2) ) * sin(2*theta[i_point])/2
                              - muxz11 * sin(4*theta[i_point])/2;
                c35_2[index] =  ( B2 * pow(cos(theta[i_point]),2) + C2 * pow(sin(theta[i_point]),2)
                                - A2 * pow(cos(theta[i_point]),2) - B2 * pow(sin(theta[i_point]),2) ) * sin(2*theta[i_point])/2
                              + muxz11 * sin(4*theta[i_point])/2;
                c55_1[index] =  ( C2 - 2*B2 + A2) * pow(sin(2*theta[i_point]),2)/4 + muxz11 * pow(cos(2*theta[i_point]),2);

            }
        }

        /*-- For v_x: rho_x , center of the volum is 01 --*/
        cal_dip_area(xvec1, zvec_h, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > nghost_z1 && iz[i_point] < nz-nghost_z2-1
              && ix[i_point] > nghost_x1 && ix[i_point] < nx-nghost_x2-1 && loc_type[i_point] != 0 ) {
                index = iz[i_point]*nx+ix[i_point];
                rho_x[index] = cal_ari_averaging(loc_type[i_point], rho10[index-nx], rho10[index-nx+1], rho10[index],
                                                 rho10[index+1], area1[i_point], area2[i_point]);
            };
        }

        /*-- For v_z: rho_z , center of the volum is 10 --*/
        cal_dip_area(xvec_h, zvec1, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer ,ix, iz, loc_type, area1, area2, theta);
        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > nghost_z1 && iz[i_point] < nz-nghost_z2-1
              && ix[i_point] > nghost_x1 && ix[i_point] < nx-nghost_x2-1 && loc_type[i_point] != 0 ) {
                index = iz[i_point]*nx+ix[i_point];
                rho_z[index] = cal_ari_averaging(loc_type[i_point], rho01[index-1], rho01[index], rho01[index+nx-1],
                                                 rho01[index+nx], area1[i_point], area2[i_point]);
            }
        }

    }


    free(ix);
    free(iz);
    free(loc_type);
    free(area1);
    free(area2);
    free(theta);
    return 0;
}


/* Calculate volume arithmetic averaging */
float cal_ari_averaging(int loc_type, float val00, float val01, float val10, float val11, float area1, float area2)
{
    float para;
    switch(loc_type) {
        case 1:
            para =  val00                 * area1 + (val11+val01+val10)/3 * area2 ;
            break;
        case 2:
            para =  (val00+val01)/2       * area1 + (val10+val11)/2       * area2 ;
            break;
        case 3:
            para =  (val00+val01+val11)/3 * area1 + val10                 * area2 ;
            break;
        case 4:
            para =  val01                 * area1 + (val00+val10+val11)/3 * area2 ;
            break;
        case 5:
            para =  (val01 + val11)/2     * area1 + (val00+val10)/2       * area2 ;
            break;
        case 6:
            para =  val11                 * area1 + (val00+val10+val01)/3 * area2 ;
            break;
        case 7:
            para =  val00                 * area2 + (val11+val01+val10)/3 * area1 ;
            break;
        case 8:
            para =  (val00+val01)/2       * area2 + (val10+val11)/2       * area1 ;
            break;
        case 9:
            para =  (val00+val01+val11)/3 * area2 + val10                 * area1 ;
            break;
        case 10:
            para =  val01                 * area2 + (val00+val10+val11)/3 * area1 ;
            break;
        case 11:
            para =  (val01+val11)/2       * area2 + (val00+val10)/2       * area1 ;
            break;
        case 12:
            para =  val11                 * area2 + (val00+val10+val01)/3 * area1 ;
            break;
        default:
            printf("There is no loc type! Please check out!");
            return;
    }

    return para/(area1+area2);
}

/* Calculate volume harmonic averaging */
float cal_har_averaging(int loc_type, float val00, float val01, float val10, float val11, float area1, float area2)
{
    float para;
    switch(loc_type) {
        case 1:
            para =  1/val00               * area1 + 3/(val11+val01+val10) * area2 ;
            break;
        case 2:
            para =  2/(val00+val01)       * area1 + 2/(val10+val11)       * area2 ;
            break;
        case 3:
            para =  3/(val00+val01+val11) * area1 + 1/val10               * area2 ;
            break;
        case 4:
            para =  1/val01               * area1 + 3/(val00+val10+val11) * area2 ;
            break;
        case 5:
            para =  2/(val01 + val11)     * area1 + 2/(val00+val10)       * area2 ;
            break;
        case 6:
            para =  1/val11               * area1 + 3/(val00+val10+val01) * area2 ;
            break;
        case 7:
            para =  1/val00               * area2 + 3/(val11+val01+val10) * area1 ;
            break;
        case 8:
            para =  2/(val00+val01)       * area2 + 2/(val10+val11)       * area1 ;
            break;
        case 9:
            para =  3/(val00+val01+val11) * area2 + 1/val10               * area1 ;
            break;
        case 10:
            para =  1/val01               * area2 + 3/(val00+val10+val11) * area1 ;
            break;
        case 11:
            para =  2/(val01+val11)/2     * area2 + 2/(val00+val10)       * area1 ;
            break;
        case 12:
            para =  1/val11               * area2 + 3/(val00+val10+val01) * area1 ;
            break;
        default:
            printf("There is no loc type! Please check out!");
            return;
    }

    return (area1+area2)/para;
}

