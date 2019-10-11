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

#include "elastic2d_src.h"
#include "elastic2d_stf.h"
#include "elastic2d_lebedev.h"
#include "staggered_fd_coef.h"
#include "elastic2d_filter.h"
#include "write_snapshots.h"
#include "write_seismograms.h"

//#include "write_seismograms.h"

int elastic2d_lebedev(float dx, float dz, int nx, int nz, int nt, float dt,
                    int half_fd_stencil, int spatial_difference_method, int filter_method,
                    float *xvec1, float *zvec1, float *xvec2, float *zvec2,
                    float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
                    float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
                    float *B01, float *B10, struct Src src, int source_impulse_method,
                    int *boundary_type, int *boundary_layer_number,
                    float *Txx_1, float *Txx_2, float *Txz_1, float *Txz_2, float *Tzz_1, float *Tzz_2,
                    float *Vx_1, float *Vx_2, float *Vz_1, float *Vz_2,
                    float *hTxx_1, float *hTxx_2, float *hTxz_1, float *hTxz_2, float *hTzz_1, float *hTzz_2,
                    float *hVx_1, float *hVx_2, float *hVz_1, float *hVz_2,
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
    ni2 = nx - half_fd_stencil-1;
    nk1 = half_fd_stencil;
    nk2 = nz - half_fd_stencil-1;
    //fprintf(stdout, "calculation area: %d %d %d %d\n", ni1, ni2, nk1, nk2);

/* ! Check stability */
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
            }
            break;
        case 2:
            for (i = 0; i < half_fd_stencil; i++) {
                fdx[i] = esgfd_coef(half_fd_stencil,i)/dx;
                fdz[i] = esgfd_coef(half_fd_stencil,i)/dz;
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
    fprintf(stdout, "RUN");
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
                            hVx_1, hVz_1, Txx_1, Txz_1, Tzz_1,
                            hVx_2, hVz_2, Txx_2, Txz_2, Tzz_2);

        /* Force source */
        ierr =  src_tti_force(hVx_1, hVz_1, hVx_2, hVz_2,
                              xvec1, zvec1,  /* for the integer grids  */
                              xvec2, zvec2,  /* for the half grids     */
                              source_impulse_method, src,
                              current_time, dt,
                              dx, dz, nx);

        /* Absorbing boundary condition */
        ierr = abs_exp_velocity(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                            boundary_type, boundary_layer_number,
                            hVx_1, hVx_2, hVz_1, hVz_2);

        /* moment equation to update velocity */
        ierr = update_velocity(nx, nz, dt, Vx_1, Vz_1, Vx_2, Vz_2,
                            hVx_1, hVz_1, hVx_2, hVz_2, B10, B01);

        /* Filtering velocity */
        if (filter_method != NOFILTER) {
            ierr = filter_velocity(filter_method, half_fd_stencil, nx,
                                ni1, ni2, nk1, nk2,
                                Vx_1, Vx_2, Vz_1, Vz_2);
        }

        /*------------------------------------------------------*/
        /* update stress (Hoke law)                             */
        /*------------------------------------------------------*/
        current_time = dt*(it+0.5);
        /* Free surface */


        /* RHS of Hook law euqation (stress-strain equation) */
        ierr = cal_hook(fdx, fdz, half_fd_stencil, ni1, ni2, nk1, nk2, nx,
                        hTxx_1, hTzz_1, hTxz_1, Vx_1 , Vz_1,
                        hTxx_2, hTzz_2, hTxz_2, Vx_2 , Vz_2,
                        c11_1 , c13_1 , c15_1 , c33_1, c35_1, c55_1,
                        c11_2 , c13_2 , c15_2 , c33_2, c35_2, c55_2 );

        /* Moment source */
        ierr = src_tti_moment(hTxx_1, hTxx_2, hTxz_1, hTxz_2, hTzz_1, hTzz_2,
                             xvec1, zvec1, xvec2, zvec2,
                             source_impulse_method, src, current_time, dt,
                             dx, dz, nx);


        /* update stress */
        ierr = update_stress(nx, nz, dt, Txx_1, Tzz_1, Txz_1, Txx_2, Tzz_2, Txz_2,
                        hTxx_1, hTzz_1, hTxz_1, hTxx_2, hTzz_2, hTxz_2);

        /* Absorbing boundary condition */
        ierr = abs_exp_stresses(ni1, ni2, nk1, nk2, nx, nz, half_fd_stencil,
                        boundary_type,  boundary_layer_number,
                        Txx_1, Txx_2, Tzz_1, Tzz_2, Txz_1, Txz_2);

        /* Filtering stresses */
        if (filter_method != NOFILTER) {
            ierr = filter_stresses(filter_method, half_fd_stencil, nx,
                            ni1, ni2, nk1, nk2,
                            Txx_1, Tzz_1, Txz_1,
                            Txx_2, Tzz_2, Txz_2);
        }


        /* receiver and snapshot */
        ierr = write_seismograms(Vx_1 , Vx_2,
                                 Vz_1 , Vz_2,
                                 Txx_1, Tzz_1,
                                 Txz_1, Txx_2,
                                 Tzz_2, Txz_2,
                                 xvec1, xvec2,
                                 zvec1, zvec2,
                                 it+1, nt, dt,
                                 nx, nz, dx, dz,
                                 nreceivers, xr, zr,
                                 save_ASCII_seismograms, save_binary_seismograms,
                                 NSTEP_BETWEEN_OUTPUT_SEISMOS, seismotype);

        ierr = write_snapshots(Vx_1,  Vx_2,
                    Vz_1,  Vz_2,
                    Txx_1, Tzz_1,
                    Txz_1, Txx_2,
                    Tzz_2, Txz_2,
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
int cal_momentum(double *fdx, double *fdz, int half_fd_stencil,
                 int ni1, int ni2, int nk1, int nk2, int nx,
                 float *hVx_1, float *hVz_1, float *Txx_1, float *Txz_1, float *Tzz_1,
                 float *hVx_2, float *hVz_2, float *Txx_2, float *Txz_2, float *Tzz_2)
{
    int ix, iz, i_diff, indx;
    float Txx_x01, Txz_z01, Txz_x01, Tzz_z01;
    float Txx_x10, Txz_z10, Txz_x10, Tzz_z10;

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx+ix;

            /*-- for the 1st set of grids --*/
            /* Vx-01 Vz-10 */
            Txx_x01 = 0.0; Txz_z01 = 0.0; Txz_x01 = 0.0; Tzz_z01 = 0.0;
            Txx_x10 = 0.0; Txz_z10 = 0.0; Txz_x10 = 0.0; Tzz_z10 = 0.0;
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {
                Txx_x01 += Dx01(Txx_1,fdx,i_diff,ix,iz,nx);
                Txz_z01 += Dz01(Txz_1,fdz,i_diff,ix,iz,nx);
                Txz_x10 += Dx10(Txz_1,fdx,i_diff,ix,iz,nx);
                Tzz_z10 += Dz10(Tzz_1,fdz,i_diff,ix,iz,nx);

                Txx_x10 += Dx10(Txx_2,fdx,i_diff,ix,iz,nx);
                Txz_z10 += Dz10(Txz_2,fdz,i_diff,ix,iz,nx);
                Txz_x01 += Dx01(Txz_2,fdx,i_diff,ix,iz,nx);
                Tzz_z01 += Dz01(Tzz_2,fdz,i_diff,ix,iz,nx);
            }

            hVx_1[indx] = Txx_x01 + Txz_z01;
            hVz_1[indx] = Txz_x10 + Tzz_z10;
            hVx_2[indx] = Txx_x10 + Txz_z10;
            hVz_2[indx] = Txz_x01 + Tzz_z01;


        }
    }

    return 0;
}


/* RHS of stress-strain equation (Hook law) */
int cal_hook(double *fdx, double *fdz, int half_fd_stencil, int ni1, int ni2, int nk1, int nk2, int nx,
             float *hTxx_1, float *hTzz_1, float *hTxz_1, float *Vx_1 , float *Vz_1,
             float *hTxx_2, float *hTzz_2, float *hTxz_2, float *Vx_2 , float *Vz_2,
             float *c11_1 , float *c13_1 , float *c15_1 , float *c33_1, float *c35_1, float *c55_1,
             float *c11_2 , float *c13_2 , float *c15_2 , float *c33_2, float *c35_2, float *c55_2)
{
    int ix, iz, i_diff, indx;
    float Vx_x00, Vz_z00, Vx_z00, Vz_x00; /* partial difference centered of (i,j) */
    float Vx_x11, Vz_z11, Vx_z11, Vz_x11; /* partial difference centered of (i+1/2,j+1/2) */

    for (iz = nk1; iz < nk2; iz++) {
        for (ix = ni1; ix < ni2; ix++) {

            indx = iz*nx + ix;

            Vx_x00 = 0.0; Vz_z00 = 0.0; Vx_z00 = 0.0; Vz_x00 = 0.0;
            Vx_x11 = 0.0; Vz_z11 = 0.0; Vx_z11 = 0.0; Vz_x11 = 0.0;
            /*-- for the 1st set of grids --*/
            for (i_diff = 1; i_diff <= half_fd_stencil; i_diff++) {

                Vx_x00 += Dx00(Vx_1,fdx,i_diff,ix,iz,nx);
                Vz_z00 += Dz00(Vz_1,fdz,i_diff,ix,iz,nx);
                Vx_z11 += Dz11(Vx_1,fdz,i_diff,ix,iz,nx);
                Vz_x11 += Dx11(Vz_1,fdx,i_diff,ix,iz,nx);

                Vx_x11 += Dx11(Vx_2,fdx,i_diff,ix,iz,nx);
                Vz_z11 += Dz11(Vz_2,fdz,i_diff,ix,iz,nx);
                Vx_z00 += Dz00(Vx_2,fdz,i_diff,ix,iz,nx);
                Vz_x00 += Dx00(Vz_2,fdx,i_diff,ix,iz,nx);

            }


            hTxx_1[indx] = c11_1[indx] * Vx_x00 + c13_1[indx] * Vz_z00 \
                         + c15_1[indx] * Vx_z00 + c15_1[indx] * Vz_x00;

            hTzz_1[indx] = c13_1[indx] * Vx_x00 + c33_1[indx] * Vz_z00 \
                         + c35_1[indx] * Vx_z00 + c35_1[indx] * Vz_x00;

            hTxz_1[indx] = c15_2[indx] * Vx_x11 + c35_2[indx] * Vz_z11 \
                         + c55_1[indx] * Vx_z11 + c55_1[indx] * Vz_x11;

            hTxx_2[indx] = c11_2[indx] * Vx_x11 + c13_2[indx] * Vz_z11 \
                         + c15_2[indx] * Vx_z11 + c15_2[indx] * Vz_x11;

            hTzz_2[indx] = c13_2[indx] * Vx_x11 + c33_2[indx] * Vz_z11 \
                         + c35_2[indx] * Vx_z11 + c35_2[indx] * Vz_x11;

            hTxz_2[indx] = c15_1[indx] * Vx_x00 + c35_1[indx] * Vz_z00 \
                         + c55_2[indx] * Vx_z00 + c55_2[indx] * Vz_x00;

        }
    }

    return 0;
}



int update_velocity(int nx, int nz, float dt,
        float *Vx_1 , float *Vz_1 , float *Vx_2 , float *Vz_2 ,
        float *hVx_1, float *hVz_1, float *hVx_2, float *hVz_2,
        float *B10, float *B01)
{
    int ix, iz, indx;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            indx = iz*nx + ix;

            Vx_1[indx] += dt*B01[indx] * hVx_1[indx] ;
            Vz_1[indx] += dt*B10[indx] * hVz_1[indx] ;

            Vx_2[indx] += dt*B10[indx] * hVx_2[indx] ;
            Vz_2[indx] += dt*B01[indx] * hVz_2[indx] ;

        }
    }

    return 0;
}

int update_stress(int nx, int nz, float dt,
                  float *Txx_1 , float *Tzz_1 , float *Txz_1,
                  float *Txx_2 , float *Tzz_2 , float *Txz_2 ,
                  float *hTxx_1, float *hTzz_1, float *hTxz_1,
                  float *hTxx_2, float *hTzz_2, float *hTxz_2 )
{
    int ix, iz, indx;
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nx; ix++) {
            indx = iz*nx + ix;

            Txx_1[indx] += dt * hTxx_1[indx];
            Tzz_1[indx] += dt * hTzz_1[indx];
            Txz_1[indx] += dt * hTxz_1[indx];

            Txx_2[indx] += dt * hTxx_2[indx];
            Tzz_2[indx] += dt * hTzz_2[indx];
            Txz_2[indx] += dt * hTxz_2[indx];

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
