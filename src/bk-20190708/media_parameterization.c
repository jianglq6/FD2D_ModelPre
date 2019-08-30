/*************************************************************************
 *
 * This function is used for media parameterization, which include many
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
 * Date: 02/2019
 *
 * History:
 *     02/2019: Original version created by Luqian Jiang
 *
 *
 ************************************************************************/

#include "global_used.h"
#include "cal_dip_area.h"
#include "media_parameterization.h"
#include "para_assign.h"

/*********************************************************************
 * For LOC, VOL, GRI:
 *   lam, lam2mu - (j, i)
 *   rho_x       - (j, i+1/2)
 *   rho_z       - (j+1/2, i)
 *   mu_xz       - (j+1/2, i+1/2)
 ********************************************************************/

/* LOC */
int media_parameterization_loc(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *xvec, float *zvec,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam, float *mu, float *muxz, float *rho_x, float *rho_z)
{
    printf("LOC: dx: %f, dz: %f, num_of_interfaces: %d\n",dx,dz,number_of_interfaces);

    int ierr, ix, iz;

    // Need error judgement !!!
    /*lam2mu, lam*/
    parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz,  MU, mu);
    parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, LAM, lam);

    /*rho_x*/
    ierr = parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin, dx, dz, nx, nz, RHO, rho_x);


    /*rho_z*/
    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin+dz/2, dx, dz, nx, nz, RHO, rho_z);


    /*muxz*/
    parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin+dz/2, dx, dz, nx, nz, MU, muxz);

    return 0;
}


/* GRI*/
int media_parameterization_gri(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *xvec, float *zvec,
        float *LAM, float *MU, float *RHO, float *L2M,
        float *lam, float *lam2mu, float *muxz,  float *Bx, float *Bz, float *rho, float *mu)
{
    printf("GRI: dx: %f, dz: %f, num_of_interfaces: %d\n",dx,dz,number_of_interfaces);

    int ix, iz;
    int ierr;


    /*lam2mu, lam*/
    parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, L2M, lam2mu);
    parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, LAM, lam);

    /* rho */
    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, RHO, rho);
    if (!ierr) {
        /* Bx */
        for (ix = 0; ix < nx-1; ix++) {
            for (iz = 0; iz < nz; iz++) {
                Bx[iz*nx + ix] = 1/rho[iz*nx + ix] + 1/rho[iz*nx + (ix+1)];
            }
        }
        for (iz = 0; iz < nz; iz++) {
            Bx[iz*nx + nx] = Bx[iz*nx + (nx-1)];
        }

        /* Bz */
        for (ix = 0; ix < nx; ix++) {
            for (iz = 1; iz < nz; iz++) {
                Bz[iz*nx + ix] = 1/rho[iz*nx + ix] + 1/rho[(iz-1)*nx + ix];
            }
        }
        for (ix = 0; ix < nx; ix++) {
            Bz[0*nx + ix] = Bz[1*nx + ix];
        }
    }

    /* muxz */
    parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, MU, mu);
    for (ix = 0; ix < nx-1; ix++) {
        for (iz = 1; iz < nz; iz++) {
            muxz[iz*nx + ix] = 4.0 / ( 1/mu[(iz-1)*nx + ix] + 1/mu[iz*nx + ix    ]
                                    +  1/mu[(iz-1)*nx + ix] + 1/mu[iz*nx + (ix+1)]);
        }
    }
    for (iz = 1; iz < nz; iz++) {
        muxz[iz*nx + nx] = mu[iz*nx + (nx-1)];
    }
    for (ix = 0; ix < nx; ix++) {
        muxz[0*nx + ix] = mu[1*nx + ix];
    }


    free(rho);
    free(mu);

    return 0;

}


int media_parameterization_vol(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *xvec, float *zvec,
        float *LAM, float *MU, float *RHO, float *lam11, float *mu00, float *mu11, float *rho01, float *rho10,
        float *lam, float *mu, float *muxz,  float *rho_x, float *rho_z)
{
    int   *ix = NULL, *iz = NULL, *loc_type = NULL;
    float *area1 = NULL, *area2 = NULL, *theta = NULL;
    float xvec_h[nx], zvec_h[nz];  // xvec-dx/2, zvec-dz/2.
    int   ix_h, iz_h, npoints_layer, i_point, index, ierr, i_interface;

    int i,j;

    for (ix_h = 0; ix_h < nx; ix_h++)
        xvec_h[ix_h] = xvec[ix_h] - dx/2;
    for (iz_h = 0; iz_h < nz; iz_h++)
        zvec_h[iz_h] = zvec[iz_h] - dz/2;

    ix = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    iz = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    loc_type = (int*)malloc(MAX_LAYER_POINTS*sizeof(int));
    area1 = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    area2 = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    theta = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));

    /* assginment first */
    // Need error judgement !!!
    /* mu, lam*/
    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, MU, mu);
    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, LAM, lam);
    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin, dx, dz, nx, nz, MU, mu00);

    /*rho_x*/
    ierr = parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin, dx, dz, nx, nz, RHO, rho_x);

    ierr = parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin, dx, dz, nx, nz, RHO, rho01);


    /*rho_z*/
    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin+dz/2, dx, dz, nx, nz, RHO, rho_z);

    ierr = parameter_assignment(interface, number_of_interfaces, xmin, zmin+dz/2, dx, dz, nx, nz, RHO, rho10);

    /*muxz*/
    ierr = parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin+dz/2, dx, dz, nx, nz, MU, muxz);

    ierr = parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin+dz/2, dx, dz, nx, nz, MU, mu11);
    ierr = parameter_assignment(interface, number_of_interfaces, xmin+dx/2, zmin+dz/2, dx, dz, nx, nz, LAM, lam11);

    /* volumn arithmetic and harmonic averaging */
    /* calculate the area */
    for (i_interface = 0; i_interface < number_of_interfaces; i_interface++) {

        /*-- For sigma_xx, sigma_zz: lam, mu, center of the volum is 00 --*/
        cal_dip_area(xvec_h, zvec_h, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer,ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > 0 && iz[i_point] < nz-1 && ix[i_point] > 0 && ix[i_point] < nx-1 ) {
                index = iz[i_point]*nx+ix[i_point];
                lam[index] = cal_para_vol(loc_type[i_point], lam11[index-1-nx], lam11[index-nx], lam11[index-1],
                                          lam11[index], area1[i_point], area2[i_point]);
                mu[index]  = cal_para_vol(loc_type[i_point],  mu11[index-1-nx],  mu11[index-nx],  mu11[index-1],
                                           mu11[index], area1[i_point], area2[i_point]);
                if (area1[i_point] < 0)
                    printf("00: area1: %d %d is smaller than 0! %f %f\n",ix[i_point], iz[i_point],area1[i_point],area2[i_point]);
                if (area2[i_point] < 0)
                    printf("00: area2: %d %d is smaller than 0! %f %f \n",ix[i_point],iz[i_point],area1[i_point],area2[i_point]);

            }

        }

        /*-- For sigma_xz: muxz , center of the volum is 11 --*/
        cal_dip_area(xvec, zvec, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > 0 && iz[i_point] < nz-1 && ix[i_point] > 0 && ix[i_point] < nx-1 ) {
                index = iz[i_point]*nx+ix[i_point];
                muxz[index] = cal_para_vol(loc_type[i_point], mu00[index], mu00[index+1], mu00[index+nx],
                                        mu00[index+nx+1], area1[i_point], area2[i_point]);
                if (area1[i_point] < 0)
                    printf("11: area1: %d %d is smaller than 0! %f \n",ix[i_point], iz[i_point],area1[i_point]);
                if (area2[i_point] < 0)
                    printf("11: area2: %d %d is smaller than 0! %f \n",ix[i_point], iz[i_point],area2[i_point]);
            }
        }

        /*-- For v_x: rho_x , center of the volum is 01 --*/
        cal_dip_area(xvec, zvec_h, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer, ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > 0 && iz[i_point] < nz-1 && ix[i_point] > 0 && ix[i_point] < nx-1 ) {
                index = iz[i_point]*nx+ix[i_point];
                rho_x[index] = cal_para_vol(loc_type[i_point], rho10[index-nx], rho10[index-nx+1], rho10[index],
                                          rho10[index+1], area1[i_point], area2[i_point]);
                if (area1[i_point] < 0)
                    printf("01: area1: %d %d is smaller than 0! %f \n",ix[i_point], iz[i_point],area1[i_point]);
                if (area2[i_point] < 0)
                    printf("01: area2: %d %d is smaller than 0! %f \n",ix[i_point], iz[i_point],area2[i_point]);
            };
        }

        /*-- For v_z: rho_z , center of the volum is 10 --*/
        cal_dip_area(xvec_h, zvec, dx, dz, nx, nz,
            interface[i_interface].x_loc, interface[i_interface].z_loc, interface[i_interface].npoints_interfaces,
            &npoints_layer ,ix, iz, loc_type, area1, area2, theta);

        for (i_point = 0; i_point < npoints_layer; i_point++) {
            if ( iz[i_point] > 0 && iz[i_point] < nz-1 && ix[i_point] > 0 && ix[i_point] < nx-1 ) {
                index = iz[i_point]*nx+ix[i_point];
                rho_z[index] = cal_para_vol(loc_type[i_point], rho01[index-1], rho01[index], rho01[index+nx-1],
                                          rho01[index+nx], area1[i_point], area2[i_point]);
                if (area1[i_point] < 0)
                    printf("10: area1: %d %d is smaller than 0! %f %f \n",ix[i_point], iz[i_point], area1[i_point], area2[i_point]);
                if (area2[i_point] < 0)
                    printf("10: area2: %d %d is smaller than 0! %f \n",ix[i_point], iz[i_point], area2[i_point]);
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

float cal_para_vol(int loc_type, float val00, float val01, float val10, float val11, float area1, float area2)
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
            para =  (val00+val10)/2       * area1 + (val01+val11)/2       * area2 ;
            break;
        case 5:
            para =  val01                 * area1 + (val00+val10+val11)/3 * area2 ;
            break;
        case 6:
            para =  (val00+val10+val01)/3 * area1 + val11                 * area2 ;
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
            para =  (val00+val10)/2       * area2 + (val01+val11)/2       * area1 ;
            break;
        case 11:
            para =  val01                 * area2 + (val00+val10+val11)/3 * area1 ;
            break;
        case 12:
            para =  (val00+val10+val01)/3 * area2 + val11                 * area1 ;
            break;
        default:
            printf("There is no loc type! Please check out!");
            return;
    }

    return para/(area1+area2);
}
//int media_parameterization_tti()
