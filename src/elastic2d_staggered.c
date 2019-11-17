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

#include "elastic2d_src.h"
#include "elastic2d_stf.h"
#include "elastic2d_staggered.h"
#include "staggered_fd_coef.h"
#include "write_snapshots.h"
#include "write_seismograms.h"

int elastic2d_staggered(float dx, float dz, int nx, int nz, int nt, float dt,
                    int half_fd_stencil, int spatial_difference_method,
                    float *xvec1, float *zvec1, float *xvec2, float *zvec2,
                    float *c11, float *c13, float *c33, float *c55,
                    float *B01, float *B10,
                    struct Src src, int source_impulse_method,
                    int *boundary_type, int *boundary_layer_number,
                    float *Txx, float *Txz, float *Tzz,
                    float *Vx, float *Vz,
                    float *hTxx, float *hTxz, float *hTzz,
                    float *hVx, float *hVz,
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

    ni1 = half_fd_stencil;
    ni2 = nx - half_fd_stencil;
    nk1 = half_fd_stencil;
    nk2 = nz - half_fd_stencil;
    printf("stencil:%d, ni1: %d, ni2: %d, nk1: %d, nk2: %d. \n",
      half_fd_stencil, ni1, ni2, nk1, nk2);

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
                            hVx, hVz, Txx, Txz, Tzz);

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
        ierr = abs_exp_velocity_ssg(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                            boundary_type, boundary_layer_number,
                            Vx, Vz);

        /*------------------------------------------------------*/
        /* update stress (Hoke law)                             */
        /*------------------------------------------------------*/
        current_time = dt*(it+0.5);
        /* Free surface */

        /* RHS of Hook law euqation (stress-strain equation) */
        ierr = cal_hook_ssg(fdx, fdz, half_fd_stencil, ni1, ni2, nk1, nk2, nx,
                        hTxx, hTzz, hTxz, Vx , Vz,
                        c11, c13, c33, c55);

        /* Moment source */
        ierr = src_moment_ssg(hTxx, hTxz, hTzz,
                             xvec1, zvec1, xvec2, zvec2,
                             source_impulse_method, src, current_time, dt,
                             dx, dz, nx);

        /* update stress */
        ierr = update_stress_ssg(nx, nz, dt, Txx, Tzz, Txz, hTxx, hTzz, hTxz);

        /* Absorbing boundary condition */
        ierr = abs_exp_stresses_ssg(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                        boundary_type, boundary_layer_number,
                        Txx, Tzz, Txz);

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


    return 0;
}


/* RHS of equation of motion */
int cal_momentum_ssg(double *fdx, double *fdz, int half_fd_stencil,
                 int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx, float *hVz, float *Txx, float *Txz, float *Tzz)
{
    int ix, iz, i_diff, indx;
    float Txx_x01, Txz_z01, Txz_x10, Tzz_z10;

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx+ix;

            /*-- for the 1st set of grids --*/
            /* Vx-01 Vz-10 */
            Txx_x01 = 0.0; Txz_z01 = 0.0; Txz_x10 = 0.0; Tzz_z10 = 0.0;
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {
                Txx_x01 += Dx01(Txx,fdx,i_diff,ix,iz,nx);
                Txz_z01 += Dz01(Txz,fdz,i_diff,ix,iz,nx);
                Txz_x10 += Dx10(Txz,fdx,i_diff,ix,iz,nx);
                Tzz_z10 += Dz10(Tzz,fdz,i_diff,ix,iz,nx);
            }

            hVx[indx] = Txx_x01 + Txz_z01;
            hVz[indx] = Txz_x10 + Tzz_z10;

        }
    }

    return 0;
}

/* RHS of stress-strain equation (Hook law) */
int cal_hook_ssg(double *fdx, double *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx, float *hTzz, float *hTxz, float *Vx , float *Vz,
             float *c11 , float *c13 , float *c33 , float *c55)
{
    int ix, iz, i_diff, indx;
    float Vx_x00, Vz_z00; /* partial difference centered of (i,j) */
    float Vx_z11, Vz_x11; /* partial difference centered of (i+1/2,j+1/2) */

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx + ix;

            Vx_x00 = 0.0; Vz_z00 = 0.0;
            Vx_z11 = 0.0; Vz_x11 = 0.0;
            /*-- for the 1st set of grids --*/
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {

                Vx_x00 += Dx00(Vx, fdx, i_diff, ix, iz, nx);
                Vz_z00 += Dz00(Vz, fdz, i_diff, ix, iz, nx);
                Vx_z11 += Dz11(Vx, fdz, i_diff, ix, iz, nx);
                Vz_x11 += Dx11(Vz, fdx, i_diff, ix, iz, nx);

            }

            hTxx[indx] = c11[indx] * Vx_x00 + c13[indx] * Vz_z00;
            hTzz[indx] = c13[indx] * Vx_x00 + c33[indx] * Vz_z00;
            hTxz[indx] = c55[indx] * Vx_z11 + c55[indx] * Vz_x11;

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
