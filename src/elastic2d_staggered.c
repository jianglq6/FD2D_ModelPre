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
 * History: 11/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include "utility_vector.h"
#include "elastic2d_src.h"
#include "elastic2d_stf.h"
#include "elastic2d_staggered.h"
#include "staggered_fd_coef.h"
#include "elastic2d_abs_exp.h"
#include "elastic2d_abs_npml.h"
#include "write_snapshots.h"
#include "write_seismograms.h"

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
                    bool use_binary_for_wavefield_dumps)
{
    int ni1, ni2, nk1, nk2, it, i, ierr, is;
    double *fdx = NULL, *fdz = NULL;
    float current_time;
    float *Ax_regu = NULL, *Bx_regu = NULL, *Dx_regu = NULL;
    float *Az_regu = NULL, *Bz_regu = NULL, *Dz_regu = NULL;
    float *Ax_half = NULL, *Bx_half = NULL, *Dx_half = NULL;
    float *Az_half = NULL, *Bz_half = NULL, *Dz_half = NULL;
    float *Txx_x1 = NULL, *Txz_x1 = NULL, *Vx_x1 = NULL, *Vz_x1 = NULL;  //x1-left
    float *Txx_x2 = NULL, *Txz_x2 = NULL, *Vx_x2 = NULL, *Vz_x2 = NULL;  //x2-right
    float *Tzz_z1 = NULL, *Txz_z1 = NULL, *Vx_z1 = NULL, *Vz_z1 = NULL;  //z1-top
    float *Tzz_z2 = NULL, *Txz_z2 = NULL, *Vx_z2 = NULL, *Vz_z2 = NULL;  //z2-bottom
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
            Txx_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
            Txz_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
             Vx_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
             Vz_x1 = creat_float_array(boundary_layer_number[0]*nz, 0.0);
        }

        if (boundary_layer_number[1] > 0) {
            Txx_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
            Txz_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
             Vx_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
             Vz_x2 = creat_float_array(boundary_layer_number[1]*nz, 0.0);
        }

        if (boundary_layer_number[2] > 0) {
            Tzz_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
            Txz_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
             Vx_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
             Vz_z1 = creat_float_array(boundary_layer_number[2]*nx, 0.0);
        }

        if (boundary_layer_number[3] > 0) {
            Tzz_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
            Txz_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
             Vx_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
             Vz_z2 = creat_float_array(boundary_layer_number[3]*nx, 0.0);
        }

        ierr = abs_npml_init(boundary_layer_number,
                         src.stf_freqfactor[0], c11, c55, B01,
                         half_fd_stencil,
                         Ax_regu, Bx_regu, Dx_regu,
                         Ax_half, Bx_half, Dx_half,
                         Az_regu, Bz_regu, Dz_regu,
                         Az_half, Bz_half, Dz_half,
                         xvec1, xvec2, zvec1, zvec2,
                         ni1, ni2, nk1, nk2, nx, dx, dz);
    }




    // TODO: Check Stability !

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

    /*==========================================================*/
    /* Time loop to simulate wavefield propagation              */
    /*==========================================================*/
    fprintf(stdout, "RUN SSG...\n");

    for (it = 0; it < nt; it++) {

        current_time = dt*it;

        if ((it+1)%NSTEP_BETWEEN_OUTPUT_INFO == 0)
            fprintf(stdout, "\n  it = %d, current time = %f \n", it+1, current_time);

        /*------------------------------------------------------*/
        /* update velocity (equation of the motion)             */
        /*------------------------------------------------------*/

        /* Free surface */

        /* RHS of equation of the motion */
        ierr = cal_momentum_ssg(fdx, fdz, half_fd_stencil,
                            ni1, ni2, nk1, nk2, nx,
                            hVx, hVz, Txx, Txz, Tzz,
                            DxTxx01, DzTxz01, DxTxz10, DzTzz10);

        /* Absorbing boundary condition */
        if (abs_type == ABS_NPML) {
            ierr = abs_npml_velocity_staggered(ni1, ni2, nk1, nk2, nx, nz,
                           boundary_layer_number, dt,
                           hVx, hVz,
                           Ax_regu, Bx_regu, Dx_regu,
                           Ax_half, Bx_half, Dx_half,
                           Az_regu, Bz_regu, Dz_regu,
                           Az_half, Bz_half, Dz_half,
                           Txx_x1, Txx_x2,
                           Tzz_z1, Tzz_z2,
                           Txz_x1, Txz_x2,
                           Txz_z1, Txz_z2,
                           DxTxx01, DzTxz01, DxTxz10, DzTzz10);
        }

        /* Force source */
        ierr =  src_force_ssg(hVx, hVz,
                              xvec1, zvec1,  /* for the integer grids  */
                              xvec2, zvec2,  /* for the half grids     */
                              source_impulse_method, src,
                              current_time, dt,
                              dx, dz, nx);

        /* moment equation to update velocity */
        ierr = update_velocity_ssg(nx, nz, dt, Vx, Vz, hVx, hVz, B10, B01);

        /* Absorbing boundary condition */
        if (abs_type == ABS_EXPONENT) {
            ierr = abs_exp_velocity_ssg(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                            boundary_layer_number,
                            Vx, Vz);
        }

        /*------------------------------------------------------*/
        /* update stress (Hoke law)                             */
        /*------------------------------------------------------*/
        current_time = dt*(it+0.5);
        /* Free surface */

        /* RHS of Hook law euqation (stress-strain equation) */
        ierr = cal_hook_ssg(fdx, fdz, half_fd_stencil, ni1, ni2, nk1, nk2, nx,
                        hTxx, hTzz, hTxz, Vx , Vz,
                        c11, c13, c33, c55,
                        DxVx00, DzVz00, DxVz11, DzVx11);

       if (abs_type == ABS_NPML) {
            ierr = abs_npml_stress_staggered(ni1, ni2, nk1, nk2, nx, nz,
                        boundary_layer_number, dt,
                        hTxx, hTxz, hTzz,
                        DxVx00, DzVz00, DzVx11, DxVz11,
                        c11, c13, c33, c55,
                        Ax_regu, Bx_regu, Dx_regu,
                        Ax_half, Bx_half, Dx_half,
                        Az_regu, Bz_regu, Dz_regu,
                        Az_half, Bz_half, Dz_half,
                        Vx_x1, Vx_x2, Vx_z1, Vx_z2,
                        Vz_x1, Vz_x2, Vz_z1, Vz_z2);
        }

        /* Moment source */
        ierr = src_moment_ssg(hTxx, hTxz, hTzz,
                             xvec1, zvec1, xvec2, zvec2,
                             source_impulse_method, src, current_time, dt,
                             dx, dz, nx);



        /* update stress */
        ierr = update_stress_ssg(nx, nz, dt, Txx, Tzz, Txz, hTxx, hTzz, hTxz);

        /* Absorbing boundary condition */
        if (abs_type == ABS_EXPONENT) {
            ierr = abs_exp_stresses_ssg(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                        boundary_layer_number,
                        Txx, Tzz, Txz);
        }

        /* receiver and snapshot */
        ierr = write_seismograms_ssg(Vx, Vz,
                                     Txx, Tzz, Txz,
                                     xvec1, xvec2,
                                     zvec1, zvec2,
                                     it+1, nt, dt,
                                     nx, nz, dx, dz,
                                     nreceivers, xr, zr,
                                     save_ASCII_seismograms, save_binary_seismograms,
                                     NSTEP_BETWEEN_OUTPUT_SEISMOS, seismotype);

        ierr = write_snapshots_ssg(Vx, Vz,
                    Txx, Tzz, Txz,
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

    fprintf(stdout, "DONE \n");

    free(fdx); free(fdz);
    free_float_array(Ax_regu); free_float_array(Az_regu);
    free_float_array(Dx_regu); free_float_array(Dz_regu);
    free_float_array(Bx_regu); free_float_array(Bz_regu);
    free_float_array(Ax_half); free_float_array(Az_half);
    free_float_array(Dx_half); free_float_array(Dz_half);
    free_float_array(Bx_half); free_float_array(Bz_half);
    free_float_array(Txx_x1); free_float_array(Txz_x1);
    free_float_array(Txx_x2); free_float_array(Txz_x2);
    free_float_array(Tzz_z1); free_float_array(Txz_z1);
    free_float_array(Tzz_z2); free_float_array(Txz_z2);
    free_float_array(Vx_x1); free_float_array(Vz_x1);
    free_float_array(Vx_x2); free_float_array(Vz_x2);
    free_float_array(Vx_z1); free_float_array(Vz_z1);
    free_float_array(Vx_z2); free_float_array(Vz_z2);

    return 0;
}


/* RHS of equation of motion */
int cal_momentum_ssg(double *fdx, double *fdz, int half_fd_stencil,
                 int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx, float *hVz, float *Txx, float *Txz, float *Tzz,
                 float *DxTxx01, float *DzTxz01, float *DxTxz10, float *DzTzz10)
{
    int ix, iz, i_diff, indx;
//    float DxTxx01[indx], DzTxz01[indx], DxTxz10[indx], DzTzz10[indx];

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx+ix;

            /*-- for the 1st set of grids --*/
            /* Vx-01 Vz-10 */
            DxTxx01[indx] = 0.0; DzTxz01[indx] = 0.0;
            DxTxz10[indx] = 0.0; DzTzz10[indx] = 0.0;
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {
                DxTxx01[indx] += Dx01(Txx,fdx,i_diff,ix,iz,nx);
                DzTxz01[indx] += Dz01(Txz,fdz,i_diff,ix,iz,nx);
                DxTxz10[indx] += Dx10(Txz,fdx,i_diff,ix,iz,nx);
                DzTzz10[indx] += Dz10(Tzz,fdz,i_diff,ix,iz,nx);
            }

            hVx[indx] = DxTxx01[indx] + DzTxz01[indx];
            hVz[indx] = DxTxz10[indx] + DzTzz10[indx];

        }
    }

    return 0;
}

/* RHS of stress-strain equation (Hook law) */
int cal_hook_ssg(double *fdx, double *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx, float *hTzz, float *hTxz, float *Vx , float *Vz,
             float *c11 , float *c13 , float *c33 , float *c55,
             float *DxVx00, float *DzVz00, float *DxVz11, float *DzVx11)
{
    int ix, iz, i_diff, indx;
//    float DxVx00, DzVz00; /* partial difference centered of (i,j) */
//    float DzVx11, DxVz11; /* partial difference centered of (i+1/2,j+1/2) */

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx + ix;

            DxVx00[indx] = 0.0; DzVz00[indx] = 0.0;
            DzVx11[indx] = 0.0; DxVz11[indx] = 0.0;
            /*-- for the 1st set of grids --*/
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {

                DxVx00[indx] += Dx00(Vx, fdx, i_diff, ix, iz, nx);
                DzVz00[indx] += Dz00(Vz, fdz, i_diff, ix, iz, nx);
                DzVx11[indx] += Dz11(Vx, fdz, i_diff, ix, iz, nx);
                DxVz11[indx] += Dx11(Vz, fdx, i_diff, ix, iz, nx);

            }

            hTxx[indx] = c11[indx] * DxVx00[indx] + c13[indx] * DzVz00[indx];
            hTzz[indx] = c13[indx] * DxVx00[indx] + c33[indx] * DzVz00[indx];
            hTxz[indx] = c55[indx] * DzVx11[indx] + c55[indx] * DxVz11[indx];

        }
    }

    return 0;
}

int update_velocity_ssg(int nx, int nz, float dt,
        float *Vx, float *Vz,
        float *hVx, float *hVz,
        float *B10, float *B01)
{
    int ix, iz, indx;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            indx = iz*nx + ix;

            Vx[indx] += dt*B01[indx] * hVx[indx] ;
            Vz[indx] += dt*B10[indx] * hVz[indx] ;

        }
    }

    return 0;
}

int update_stress_ssg(int nx, int nz, float dt,
                  float *Txx , float *Tzz , float *Txz,
                  float *hTxx, float *hTzz, float *hTxz)
{
    int ix, iz, indx;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            indx = iz*nx + ix;

            Txx[indx] += dt * hTxx[indx];
            Tzz[indx] += dt * hTzz[indx];
            Txz[indx] += dt * hTxz[indx];

        }
    }

    return 0;
}
