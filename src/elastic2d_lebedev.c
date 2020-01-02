/***************************************************************************
 *
 * This function is used to simulate 2D elastic waveform in the TTI media
 *   using the high-order Lebedev scheme
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 08/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include "utility_vector.h"
#include "elastic2d_src.h"
#include "elastic2d_stf.h"
#include "elastic2d_lebedev.h"
#include "staggered_fd_coef.h"
#include "elastic2d_filter.h"
#include "elastic2d_abs_npml.h"
#include "elastic2d_abs_exp.h"
#include "write_snapshots.h"
#include "write_seismograms.h"


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
                    bool use_binary_for_wavefield_dumps)
{

    int ni1, ni2, nk1, nk2, it, i, ierr, is;
    double *fdx = NULL, *fdz = NULL;
    float current_time;
    float *Ax_regu = NULL, *Bx_regu = NULL, *Dx_regu = NULL;
    float *Az_regu = NULL, *Bz_regu = NULL, *Dz_regu = NULL;
    float *Ax_half = NULL, *Bx_half = NULL, *Dx_half = NULL;
    float *Az_half = NULL, *Bz_half = NULL, *Dz_half = NULL;
    float *Txx01_x1 = NULL, *Txz01_x1 = NULL, *Vx00_x1 = NULL, *Vz00_x1 = NULL;  //x1-left
    float *Txx01_x2 = NULL, *Txz01_x2 = NULL, *Vx00_x2 = NULL, *Vz00_x2 = NULL;  //x2-right
    float *Tzz01_z1 = NULL, *Txz01_z1 = NULL, *Vx00_z1 = NULL, *Vz00_z1 = NULL;  //z1-top
    float *Tzz01_z2 = NULL, *Txz01_z2 = NULL, *Vx00_z2 = NULL, *Vz00_z2 = NULL;  //z2-bottom
    float *Txx10_x1 = NULL, *Txz10_x1 = NULL, *Vx11_x1 = NULL, *Vz11_x1 = NULL;  //x1-left
    float *Txx10_x2 = NULL, *Txz10_x2 = NULL, *Vx11_x2 = NULL, *Vz11_x2 = NULL;  //x2-right
    float *Tzz10_z1 = NULL, *Txz10_z1 = NULL, *Vx11_z1 = NULL, *Vz11_z1 = NULL;  //z1-top
    float *Tzz10_z2 = NULL, *Txz10_z2 = NULL, *Vx11_z2 = NULL, *Vz11_z2 = NULL;  //z2-bottom
    float *fc_pml;

    /* calculation area */
    ni1 = half_fd_stencil;
    ni2 = nx - half_fd_stencil;
    nk1 = half_fd_stencil;
    nk2 = nz - half_fd_stencil;

    /* prepare for ade cfs-pml */
    if (abs_type == ABS_NPML) {

        Ax_regu = creat_float_array(nx, 0.0);
        Bx_regu = creat_float_array(nx, 1.0);
        Dx_regu = creat_float_array(nx, 0.0);
        Ax_half = creat_float_array(nx, 0.0);
        Bx_half = creat_float_array(nx, 1.0);
        Dx_half = creat_float_array(nx, 0.0);

        Az_regu = creat_float_array(nz, 0.0);
        Bz_regu = creat_float_array(nz, 1.0);
        Dz_regu = creat_float_array(nz, 0.0);
        Az_half = creat_float_array(nz, 0.0);
        Bz_half = creat_float_array(nz, 1.0);
        Dz_half = creat_float_array(nz, 0.0);

        if (boundary_layer_number[0] > 0) {
            Txx01_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
            Txz01_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
             Vx11_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
             Vz11_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
            Txx10_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
            Txz10_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
             Vx00_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
             Vz00_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
        }

        if (boundary_layer_number[1] > 0) {
            Txx01_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
            Txz01_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
             Vx11_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
             Vz11_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
            Txx10_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
            Txz10_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
             Vx00_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
             Vz00_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
        }

        if (boundary_layer_number[2] > 0) {
            Tzz01_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
            Txz01_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
             Vx11_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
             Vz11_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
            Tzz10_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
            Txz10_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
             Vx00_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
             Vz00_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
        }

        if (boundary_layer_number[3] > 0) {
            Tzz01_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
            Txz01_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
             Vx11_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
             Vz11_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
            Tzz10_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
            Txz10_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
             Vx00_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
             Vz00_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
        }

        ierr = abs_npml_init(boundary_layer_number,
                             src.stf_freqfactor[0], c11_1, c55_1, B01,
                             half_fd_stencil,
                             Ax_regu, Bx_regu, Dx_regu,
                             Ax_half, Bx_half, Dx_half,
                             Az_regu, Bz_regu, Dz_regu,
                             Az_half, Bz_half, Dz_half,
                             xvec1, xvec2, zvec1, zvec2,
                             ni1, ni2, nk1, nk2, nx, dx, dz);

    }


    //fprintf(stdout, "calculation area: %d %d %d %d\n", ni1, ni2, nk1, nk2);


/* ! TODO: Check stability */
    //

    /*---------------------------------------------------------*/
    /* staggered grid differential coefficient                 */
    /*---------------------------------------------------------*/
    fdx = (double*) malloc(half_fd_stencil*sizeof(double));
    fdz = (double*) malloc(half_fd_stencil*sizeof(double));

    switch(spatial_difference_method) {
        case 1:
            for (i = 0; i < half_fd_stencil; i++) {
                fdx[i] = ssg_coef(half_fd_stencil,i)/dx;
                fdz[i] = ssg_coef(half_fd_stencil,i)/dz;
                printf("%f %f\n",fdx[i]*dx, fdz[i]*dz );

            }
            break;
        case 2:
            for (i = 0; i < half_fd_stencil; i++) {
                fdx[i] = esgfd_coef(half_fd_stencil,i)/dx;
                fdz[i] = esgfd_coef(half_fd_stencil,i)/dz;
                printf("%f %f\n",fdx[i], fdz[i] );
            }
            break;

    }

    /* prepare stf */
    //printf("dt: %f\n",dt );
    //for (is = 0; is < src.number_of_src; is++) {
    //    ierr = prepare_source_time_function(nt, dt, src.stf_timefactor[is],
    //            src.stf_freqfactor[is], src.stf_type_id[is], is+1);
    //}

    /* TODO! check stability */

    /* Initial */

    /*==========================================================*/
    /* Time loop to simulate wavefield propagation              */
    /*==========================================================*/
    fprintf(stdout, "RUN...\n");

    for (it = 0; it < nt; it++) {

        current_time = dt*it;

        if ((it+1)%NSTEP_BETWEEN_OUTPUT_INFO == 0)
            fprintf(stdout, "\n  it = %d, current time = %f \n", it+1, current_time);

        /*------------------------------------------------------*/
        /* update velocity (equation of the motion)             */
        /*------------------------------------------------------*/

        /* Free surface */

        /* RHS of equation of the motion */
        ierr = cal_momentum(fdx, fdz, half_fd_stencil,
                            ni1, ni2, nk1, nk2, nx,
                            hVx01, hVz10, Txx00, Txz11, Tzz00,
                            hVx10, hVz01, Txx11, Txz00, Tzz11,
                            DxTxx01, DzTxz01, DxTxz01, DzTzz01,
                            DxTxx10, DzTxz10, DxTxz10, DzTzz10);

        /* Force source */
        ierr =  src_force_lebedev(hVx01, hVz10, hVx10, hVz01,
                              xvec1, zvec1,  /* for the integer grids  */
                              xvec2, zvec2,  /* for the half grids     */
                              source_impulse_method, src,
                              current_time, dt,
                              dx, dz, nx);

        /* Absorbing boundary condition */
        if (abs_type == ABS_NPML) {
            ierr = abs_npml_velocity_tti(ni1, ni2, nk1, nk2, nx, nz,
                           boundary_layer_number, dt,
                           hVx01, hVx10, hVz10, hVz01,
                           Ax_regu, Bx_regu, Dx_regu,
                           Ax_half, Bx_half, Dx_half,
                           Az_regu, Bz_regu, Dz_regu,
                           Az_half, Bz_half, Dz_half,
                           Txx01_x1, Txx01_x2,
                           Tzz01_z1, Tzz01_z2,
                           Txz01_x1, Txz01_x2,
                           Txz01_z1, Txz01_z2,
                           Txx10_x1, Txx10_x2,
                           Tzz10_z1, Tzz10_z2,
                           Txz10_x1, Txz10_x2,
                           Txz10_z1, Txz10_z2,
                           DxTxx01, DzTxz01, DxTxz01, DzTzz01,
                           DxTxx10, DzTxz10, DxTxz10, DzTzz10);
        }

        /* moment equation to update velocity */
        ierr = update_velocity(nx, nz, dt, Vx01, Vz10, Vx10, Vz01,
                            hVx01, hVz10, hVx10, hVz01, B10, B01);

        /* Absorbing boundary condition */
        if (abs_type == ABS_EXPONENT) {
            ierr = abs_exp_velocity(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                            boundary_layer_number,
                            Vx01, Vx10, Vz10, Vz01);
        }

        /* Filtering velocity */
        if (filter_method != NOFILTER) {
            ierr = filter_velocity(filter_method, half_fd_stencil, nx,
                                ni1, ni2, nk1, nk2,
                                Vx01, Vx10, Vz10, Vz01);
        }

        /*------------------------------------------------------*/
        /* update stress (Hoke law)                             */
        /*------------------------------------------------------*/
        current_time = dt*(it+0.5);
        /* Free surface */


        /* RHS of Hook law euqation (stress-strain equation) */
        ierr = cal_hook(fdx, fdz, half_fd_stencil, ni1, ni2, nk1, nk2, nx,
                        hTxx00, hTzz00, hTxz11, Vx01 , Vz10,
                        hTxx11, hTzz11, hTxz00, Vx10 , Vz01,
                        DxVx00, DzVz00, DzVx00, DxVz00,
                        DxVx11, DzVz11, DzVx11, DxVz11,
                        c11_1 , c13_1 , c15_1 , c33_1, c35_1, c55_1,
                        c11_2 , c13_2 , c15_2 , c33_2, c35_2, c55_2);

        if (abs_type == ABS_NPML) {
            ierr = abs_npml_stress_tti(ni1, ni2, nk1, nk2, nx, nz,
                        boundary_layer_number, dt,
                        hTxx00, hTxx11,
                        hTxz00, hTxz11,
                        hTzz00, hTzz11,
                        DxVx00, DzVz00, DzVx00, DxVz00,
                        DxVx11, DzVz11, DzVx11, DxVz11,
                        c11_1 , c13_1 , c15_1 , c33_1, c35_1, c55_1,
                        c11_2 , c13_2 , c15_2 , c33_2, c35_2, c55_2,
                        Ax_regu, Bx_regu, Dx_regu,
                        Ax_half, Bx_half, Dx_half,
                        Az_regu, Bz_regu, Dz_regu,
                        Az_half, Bz_half, Dz_half,
                        Vx00_x1, Vx00_x2, Vx11_x1, Vx11_x2,
                        Vx00_z1, Vx00_z2, Vx11_z1, Vx11_z2,
                        Vz00_x1, Vz00_x2, Vz11_x1, Vz11_x2,
                        Vz00_z1, Vz00_z2, Vz11_z1, Vz11_z2);
        }

        /* Moment source */
        ierr = src_moment_lebedev(hTxx00, hTxx11, hTxz11, hTxz00, hTzz00, hTzz11,
                             xvec1, zvec1, xvec2, zvec2,
                             source_impulse_method, src, current_time, dt,
                             dx, dz, nx);

        /* update stress */
        ierr = update_stress(nx, nz, dt, Txx00, Tzz00, Txz11, Txx11, Tzz11, Txz00,
                        hTxx00, hTzz00, hTxz11, hTxx11, hTzz11, hTxz00);

        /* Absorbing boundary condition */
        if (abs_type == ABS_EXPONENT) {
            ierr = abs_exp_stresses(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                        boundary_layer_number,
                        Txx00, Txx11, Tzz00, Tzz11, Txz11, Txz00);
        }

        /* Filtering stresses */
        if (filter_method != NOFILTER) {
            ierr = filter_stresses(filter_method, half_fd_stencil, nx,
                            ni1, ni2, nk1, nk2,
                            Txx00, Tzz00, Txz11,
                            Txx11, Tzz11, Txz00);
        }


        /* receiver and snapshot */
        ierr = write_seismograms(Vx01 , Vx10,
                                 Vz10 , Vz01,
                                 Txx00, Tzz00,
                                 Txz11, Txx11,
                                 Tzz11, Txz00,
                                 xvec1, xvec2,
                                 zvec1, zvec2,
                                 it+1, nt, dt,
                                 nx, nz, dx, dz,
                                 nreceivers, xr, zr,
                                 save_ASCII_seismograms, save_binary_seismograms,
                                 NSTEP_BETWEEN_OUTPUT_SEISMOS, seismotype);

        ierr = write_snapshots(Vx01,  Vx10,
                    Vz10,  Vz01,
                    Txx00, Tzz00,
                    Txz11, Txx11,
                    Tzz11, Txz00,
                    xvec1, xvec2,
                    zvec1, zvec2,
                    it+1, nt, dt,   // it = it+1
                    nx, nz, dx, dz, half_fd_stencil,
                    boundary_layer_number,
                    NSTEP_BETWEEN_OUTPUT_IMAGES,
                    output_postscript_snapshot,
                    imagetype_postscript, meshvect, modelvect, boundvect,
                    sizemax_arrows, US_LETTER,
                    output_wavefield_dumps, imagetype_wavefield_dumps,
                    use_binary_for_wavefield_dumps);


    }

    free(fdx); free(fdz);
    //if (abs_type == ABS_NPML) {
        free_float_array(Ax_regu); free_float_array(Az_regu);
        free_float_array(Dx_regu); free_float_array(Dz_regu);
        free_float_array(Bx_regu); free_float_array(Bz_regu);
        free_float_array(Ax_half); free_float_array(Az_half);
        free_float_array(Dx_half); free_float_array(Dz_half);
        free_float_array(Bx_half); free_float_array(Bz_half);
        free_float_array(Txx10_x1); free_float_array(Txz10_x1);
        free_float_array(Txx10_x2); free_float_array(Txz10_x2);
        free_float_array(Tzz10_z1); free_float_array(Txz10_z1);
        free_float_array(Tzz10_z2); free_float_array(Txz10_z2);
        free_float_array(Vx11_x1); free_float_array(Vz11_x1);
        free_float_array(Vx11_x2); free_float_array(Vz11_x2);
        free_float_array(Vx11_z1); free_float_array(Vz11_z1);
        free_float_array(Vx11_z2); free_float_array(Vz11_z2);
        free_float_array(Txx01_x1); free_float_array(Txz01_x1);
        free_float_array(Txx01_x2); free_float_array(Txz01_x2);
        free_float_array(Tzz01_z1); free_float_array(Txz01_z1);
        free_float_array(Tzz01_z2); free_float_array(Txz01_z2);
        free_float_array(Vx00_x1); free_float_array(Vz00_x1);
        free_float_array(Vx00_x2); free_float_array(Vz00_x2);
        free_float_array(Vx00_z1); free_float_array(Vz00_z1);
        free_float_array(Vx00_z2); free_float_array(Vz00_z2);
    //}

    fprintf(stdout, "DONE \n");

    return 0;
}

/* RHS of equation of motion */
int cal_momentum(double *fdx, double *fdz, int half_fd_stencil,
                 int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx01, float *hVz10, float *Txx00, float *Txz11, float *Tzz00,
                 float *hVx10, float *hVz01, float *Txx11, float *Txz00, float *Tzz11,
                 float *DxTxx01, float *DzTxz01, float *DxTxz01, float *DzTzz01,
                 float *DxTxx10, float *DzTxz10, float *DxTxz10, float *DzTzz10)
{
    int ix, iz, i_diff, indx;
//    float DxTxx01, DzTxz01, DxTxz01, DzTzz01;
//    float DxTxx10, DzTxz10, DxTxz10, DzTzz10;

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx+ix;

            /*-- for the 1st set of grids --*/
            /* Vx-01 Vz-10 */
            DxTxx01[indx] = 0.0; DzTxz01[indx] = 0.0; DxTxz01[indx] = 0.0; DzTzz01[indx] = 0.0;
            DxTxx10[indx] = 0.0; DzTxz10[indx] = 0.0; DxTxz10[indx] = 0.0; DzTzz10[indx] = 0.0;
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {
                DxTxx01[indx] += Dx01(Txx00,fdx,i_diff,ix,iz,nx);
                DzTxz01[indx] += Dz01(Txz11,fdz,i_diff,ix,iz,nx);
                DxTxz10[indx] += Dx10(Txz11,fdx,i_diff,ix,iz,nx);
                DzTzz10[indx] += Dz10(Tzz00,fdz,i_diff,ix,iz,nx);

                DxTxx10[indx] += Dx10(Txx11,fdx,i_diff,ix,iz,nx);
                DzTxz10[indx] += Dz10(Txz00,fdz,i_diff,ix,iz,nx);
                DxTxz01[indx] += Dx01(Txz00,fdx,i_diff,ix,iz,nx);
                DzTzz01[indx] += Dz01(Tzz11,fdz,i_diff,ix,iz,nx);
            }

            hVx01[indx] = DxTxx01[indx] + DzTxz01[indx];
            hVz10[indx] = DxTxz10[indx] + DzTzz10[indx];
            hVx10[indx] = DxTxx10[indx] + DzTxz10[indx];
            hVz01[indx] = DxTxz01[indx] + DzTzz01[indx];


        }
    }

    return 0;
}


/* RHS of stress-strain equation (Hook law) */
int cal_hook(double *fdx, double *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx00, float *hTzz00, float *hTxz11, float *Vx01 , float *Vz10,
             float *hTxx11, float *hTzz11, float *hTxz00, float *Vx10 , float *Vz01,
             float *DxVx00, float *DzVz00, float *DzVx00, float *DxVz00,
             float *DxVx11, float *DzVz11, float *DzVx11, float *DxVz11,
             float *c11_1 , float *c13_1 , float *c15_1 , float *c33_1, float *c35_1, float *c55_1,
             float *c11_2 , float *c13_2 , float *c15_2 , float *c33_2, float *c35_2, float *c55_2)
{
    int ix, iz, i_diff, indx;
//    float DxVx00, DzVz00, DzVx00, DxVz00; /* partial difference centered of (i,j) */
//    float DxVx11, DzVz11, DzVx11, DxVz11; /* partial difference centered of (i+1/2,j+1/2) */

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx + ix;

            DxVx00[indx] = 0.0; DzVz00[indx] = 0.0; DzVx00[indx] = 0.0; DxVz00[indx] = 0.0;
            DxVx11[indx] = 0.0; DzVz11[indx] = 0.0; DzVx11[indx] = 0.0; DxVz11[indx] = 0.0;
            /*-- for the 1st set of grids --*/
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {

                DxVx00[indx] += Dx00(Vx01,fdx,i_diff,ix,iz,nx);
                DzVz00[indx] += Dz00(Vz10,fdz,i_diff,ix,iz,nx);
                DzVx11[indx] += Dz11(Vx01,fdz,i_diff,ix,iz,nx);
                DxVz11[indx] += Dx11(Vz10,fdx,i_diff,ix,iz,nx);

                DxVx11[indx] += Dx11(Vx10,fdx,i_diff,ix,iz,nx);
                DzVz11[indx] += Dz11(Vz01,fdz,i_diff,ix,iz,nx);
                DzVx00[indx] += Dz00(Vx10,fdz,i_diff,ix,iz,nx);
                DxVz00[indx] += Dx00(Vz01,fdx,i_diff,ix,iz,nx);

            }


            hTxx00[indx] = c11_1[indx] * DxVx00[indx] + c13_1[indx] * DzVz00[indx] \
                         + c15_1[indx] * DzVx00[indx] + c15_1[indx] * DxVz00[indx];

            hTzz00[indx] = c13_1[indx] * DxVx00[indx] + c33_1[indx] * DzVz00[indx] \
                         + c35_1[indx] * DzVx00[indx] + c35_1[indx] * DxVz00[indx];

            hTxz11[indx] = c15_2[indx] * DxVx11[indx] + c35_2[indx] * DzVz11[indx] \
                         + c55_1[indx] * DzVx11[indx] + c55_1[indx] * DxVz11[indx];

            hTxx11[indx] = c11_2[indx] * DxVx11[indx] + c13_2[indx] * DzVz11[indx] \
                         + c15_2[indx] * DzVx11[indx] + c15_2[indx] * DxVz11[indx];

            hTzz11[indx] = c13_2[indx] * DxVx11[indx] + c33_2[indx] * DzVz11[indx] \
                         + c35_2[indx] * DzVx11[indx] + c35_2[indx] * DxVz11[indx];

            hTxz00[indx] = c15_1[indx] * DxVx00[indx] + c35_1[indx] * DzVz00[indx] \
                         + c55_2[indx] * DzVx00[indx] + c55_2[indx] * DxVz00[indx];

        }
    }

    return 0;
}



int update_velocity(int nx, int nz, float dt,
        float *Vx01 , float *Vz10 , float *Vx10 , float *Vz01 ,
        float *hVx01, float *hVz10, float *hVx10, float *hVz01,
        float *B10, float *B01)
{
    int ix, iz, indx;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            indx = iz*nx + ix;

            Vx01[indx] += dt*B01[indx] * hVx01[indx] ;
            Vz10[indx] += dt*B10[indx] * hVz10[indx] ;

            Vx10[indx] += dt*B10[indx] * hVx10[indx] ;
            Vz01[indx] += dt*B01[indx] * hVz01[indx] ;

        }
    }

    return 0;
}

int update_stress(int nx, int nz, float dt,
                  float *Txx00 , float *Tzz00 , float *Txz11,
                  float *Txx11 , float *Tzz11 , float *Txz00,
                  float *hTxx00, float *hTzz00, float *hTxz11,
                  float *hTxx11, float *hTzz11, float *hTxz00)
{
    int ix, iz, indx;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            indx = iz*nx + ix;

            Txx00[indx] += dt * hTxx00[indx];
            Tzz00[indx] += dt * hTzz00[indx];
            Txz11[indx] += dt * hTxz11[indx];

            Txx11[indx] += dt * hTxx11[indx];
            Tzz11[indx] += dt * hTzz11[indx];
            Txz00[indx] += dt * hTxz00[indx];

        }
    }

    return 0;
}

//float Dx00( float *a, double *fdx, int ith, int indx, int nx)
//{
//
//    float dx00;
//    dx00 = fdx[ith-1] * ( a[indx+ith-1] - a[indx-ith ] );
//    return dx00;
//}
//
//float Dx10( float *a, double *fdx, int ith, int indx, int nx)
//{
//    float dx10;
//    dx10 = fdx[ith-1] * ( a[indx+ith-1] - a[indx-ith] );
//    return dx10;
//}
//
//float Dx01( float *a, double *fdx, int ith, int indx, int nx)
//{
//    float dx01;
//    dx01 = fdx[ith-1] * ( a[indx+ith] - a[indx-ith+1] );
//    return dx01;
//}
//
//float Dx11( float *a, double *fdx, int ith, int indx, int nx)
//{
//    float dx11;
//    dx11 = fdx[ith-1] * ( a[indx+ith] - a[indx-ith+1] );
//    return dx11;
//}
//
//float Dz00( float *a, double *fdz, int ith, int indx, int nx)
//{
//    float dz00;
//    dz00 = fdz[ith-1] * ( a[indx+(ith-1)*nx] - a[indx-ith*nx] );
//    return dz00;
//}
//
//float Dz01( float *a, double *fdz, int ith, int indx, int nx)
//{
//    float dz01;
//    dz01 = fdz[ith-1] * ( a[indx+(ith-1)*nx] - a[indx-ith*nx] );
//    return dz01;
//}
//
//float Dz10( float *a, double *fdz, int ith, int indx, int nx)
//{
//    float dz10;
//    dz10 = fdz[ith-1] * ( a[indx+ith*nx] - a[indx-(ith-1)*nx] );
//    return dz10;
//}
//
//float Dz11( float *a, double *fdz, int ith, int indx, int nx)
//{
//    float dz11;
//    dz11 = fdz[ith-1] * ( a[indx+ith*nx] - a[indx-(ith-1)*nx] );
//    return dz11;
//}
//
//
