#include "elastic2d_src.h"
#include "elastic2d_lebedev.h"
#include "elastic2d_staggered.h"
#include "share_param.h"
#include "read_config_para.h"
#include "pre_model_prepare.h"

float *creat_array(int n, float v);

int main()
{
    char config_file[] = "./DATA/Par_file";
    int ierr = 0;

/*======================================================================
 *
 * read and initialize (The order is the same as the configuration file)
 *
 * ====================================================================*/

/*==================== simulation input parameters ====================*/
/*-- Time difference parameter --*/
    int nt;
    float dt;
/*-- spatial difference information --*/
    int half_fd_stencil, spatial_difference_method;
/*-- filter information (for Lebedev test) --*/
    int filter_method;

/*============ geometry and the grid setting of the model =============*/
/*-- read from Par_file --*/
    int  nx, nz, ix, iz, indx;
    float  dx, dz, xmin, zmin;
/*-- grid info after add boundary layer and FD stencil --*/
    int nx_all, nz_all;
    float xmin_all, zmin_all;
    float *xvec1, *xvec2, *zvec1, *zvec2;
    int nghost_x1, nghost_x2, nghost_z1, nghost_z2;

/*========================= source information ========================*/
    int source_impulse_method, is;
    struct Src src;      //include number of source, source location, and some other source information

/*====================== receiver information =========================*/
    int nreceivers, seismotype, NSTEP_BETWEEN_OUTPUT_SEISMOS;
    bool save_binary_seismograms, save_ASCII_seismograms;
    float *xr=NULL, *zr=NULL;

/*======================== boundary conditions ========================*/
    int boundary_type[4], boundary_layer_number[4], ib;

/*=====================  velocity and density model ===================*/
    bool use_existing_model;
    int nbmodels, save_model;
    char materialfile[MAX_VAL_LEN], interfacesfile[MAX_VAL_LEN];
    //char model_file_type[MAX_VAL_LEN], model_file[MAX_VAL_LEN];
/*============= effective media parameterization method ==============*/
    int effective_para_method;
    int prepared_media;  // Choosing SSG or Lebedev, because of media properties after model prepare.
    float *c11_1=NULL, *c13_1=NULL, *c15_1=NULL, *c33_1=NULL, *c35_1=NULL, *c55_1=NULL;
    float *c11_2=NULL, *c13_2=NULL, *c15_2=NULL, *c33_2=NULL, *c35_2=NULL, *c55_2=NULL;
    float *B01=NULL, *B10=NULL;
    // For TTI, Assignment temporary use
    float *lam2mu00=NULL, *lam2mu11=NULL, *lam11=NULL, *lam00=NULL, *mu00=NULL, *mu11=NULL;
    float *rho_x=NULL, *rho_z=NULL, *rho01=NULL, *rho10=NULL;

/*==================== display parameters ============================*/
    int NSTEP_BETWEEN_OUTPUT_INFO;

/*================ movies snapshot information =======================*/
    int NSTEP_BETWEEN_OUTPUT_IMAGES;
    bool output_wavefield_dumps, use_binary_for_wavefield_dumps;
    int imagetype_wavefield_dumps=1, imagetype_postscript=2;
    bool meshvect, modelvect, boundvect, US_LETTER, output_postscript_snapshot;
    float sizemax_arrows=1.0;

/*======================== alloc for Lebedev =========================*/
    float *Txx_1=NULL, *Txx_2=NULL, *Tzz_1=NULL, *Tzz_2=NULL, *Txz_1=NULL, *Txz_2=NULL;
    float *Vx_1=NULL,  *Vx_2=NULL,  *Vz_1=NULL,  *Vz_2=NULL;
    float *hTxx_1=NULL, *hTxx_2=NULL, *hTzz_1=NULL, *hTzz_2=NULL, *hTxz_1=NULL, *hTxz_2=NULL;
    float *hVx_1=NULL,  *hVx_2=NULL,  *hVz_1=NULL,  *hVz_2=NULL;

/*===================== get the info and print =======================*/
    ierr = get_config_info(config_file, &nt, &dt,
                           &half_fd_stencil, &spatial_difference_method,
                           &xmin, &dx, &nx, &zmin, &dz, &nz, &filter_method,
                           &source_impulse_method, &src,
                           &seismotype, &NSTEP_BETWEEN_OUTPUT_SEISMOS,
                           &save_ASCII_seismograms, &save_binary_seismograms,
                           &nreceivers, &xr, &zr,
                            boundary_type, boundary_layer_number,
                           &use_existing_model, &nbmodels,
                           materialfile, interfacesfile, &save_model,
                           &effective_para_method,
                           &NSTEP_BETWEEN_OUTPUT_INFO,
                           &NSTEP_BETWEEN_OUTPUT_IMAGES, &output_postscript_snapshot,
                           &imagetype_postscript, &meshvect, &modelvect, &boundvect,
                           &sizemax_arrows, &US_LETTER,
                           &output_wavefield_dumps, &imagetype_wavefield_dumps,
                           &use_binary_for_wavefield_dumps);


/*================== Initialize material array =======================*/
    /* prepare for boundary condition, TODO (free surface) */
    for (ib = 0; ib < 4; ib++){
        if (boundary_type[ib] == 0)
            boundary_layer_number[ib] = 0;
    }

    nghost_x1 = boundary_layer_number[0] + half_fd_stencil;
    nghost_x2 = boundary_layer_number[1] + half_fd_stencil;
    nghost_z1 = boundary_layer_number[2] + half_fd_stencil;
    nghost_z2 = boundary_layer_number[3] + half_fd_stencil;


    /* reconstruct the coordinate */
    xmin_all = xmin - (boundary_layer_number[0]+half_fd_stencil) * dx;
    zmin_all = zmin - (boundary_layer_number[2]+half_fd_stencil) * dz;
    nx_all = nx + boundary_layer_number[0] + boundary_layer_number[1] + half_fd_stencil*2;
    nz_all = nz + boundary_layer_number[2] + boundary_layer_number[3] + half_fd_stencil*2;


    printf("%d %d %d %d, nx_all: %d, nz_all: %d\n",
            nghost_x1, nghost_x2, nghost_z1, nghost_z2, nx_all, nz_all);

    zvec1 = (float*) malloc(nz_all*sizeof(float));
    zvec2 = (float*) malloc(nz_all*sizeof(float));
    xvec1 = (float*) malloc(nx_all*sizeof(float));
    xvec2 = (float*) malloc(nx_all*sizeof(float));

    /* alloc the material array */
    c11_1 = creat_array(nx_all*nz_all, -1.0);
    c13_1 = creat_array(nx_all*nz_all, -1.0);
    c15_1 = creat_array(nx_all*nz_all,  0.0);
    c33_1 = creat_array(nx_all*nz_all, -1.0);
    c35_1 = creat_array(nx_all*nz_all,  0.0);
    c55_1 = creat_array(nx_all*nz_all, -1.0);

    c11_2 = creat_array(nx_all*nz_all, -1.0);
    c13_2 = creat_array(nx_all*nz_all, -1.0);
    c15_2 = creat_array(nx_all*nz_all,  0.0);
    c33_2 = creat_array(nx_all*nz_all, -1.0);
    c35_2 = creat_array(nx_all*nz_all,  0.0);
    c55_2 = creat_array(nx_all*nz_all, -1.0);

    B10   = creat_array(nx_all*nz_all, 1.0);
    B01   = creat_array(nx_all*nz_all, 1.0);

    lam2mu00 = creat_array(nx_all*nz_all, -1.0);
    lam2mu11 = creat_array(nx_all*nz_all, -1.0);
    lam11    = creat_array(nx_all*nz_all, -1.0);
    lam00    = creat_array(nx_all*nz_all, -1.0);
    mu00     = creat_array(nx_all*nz_all, -1.0);
    mu11     = creat_array(nx_all*nz_all, -1.0);
    rho_x    = creat_array(nx_all*nz_all, -1.0);
    rho_z    = creat_array(nx_all*nz_all, -1.0);
    rho01    = creat_array(nx_all*nz_all, -1.0);
    rho10    = creat_array(nx_all*nz_all, -1.0);

    /* Grid vectors, xvec2, zvec2 are half grid */
    for (ix=0; ix<nx_all; ix++) {
        xvec1[ix] = xmin_all      + dx*ix;
        xvec2[ix] = xmin_all+dx/2 + dx*ix;
    }

    for (iz=0; iz<nz_all; iz++) {
        zvec1[iz] = zmin_all      + dz*iz;
        zvec2[iz] = zmin_all+dz/2 + dz*iz;
    }

/*=========== Model prepare: effective media parameterization =========*/
    if (!use_existing_model) {
        ierr = model_prepare(nbmodels, materialfile, interfacesfile,
                effective_para_method, save_model, &prepared_media,
                nghost_x1, nghost_x2, nghost_z1, nghost_z2,
                xvec1, zvec1, xvec2, zvec2, nx_all, nz_all, dx, dz,
                c11_1, c13_1, c15_1, c33_1, c35_1, c55_1,
                c11_2, c13_2, c15_2, c33_2, c35_2, c55_2,
                B01, B10,
                lam2mu00, lam2mu11, lam00, lam11, mu00, mu11,
                rho01, rho10, rho_x, rho_z);
    }


    /* alloc the components */
    Txx_1  = creat_array(nx_all*nz_all, 0.0);
    Txx_2  = creat_array(nx_all*nz_all, 0.0);
    Tzz_1  = creat_array(nx_all*nz_all, 0.0);
    Tzz_2  = creat_array(nx_all*nz_all, 0.0);
    Txz_1  = creat_array(nx_all*nz_all, 0.0);
    Txz_2  = creat_array(nx_all*nz_all, 0.0);
    Vz_1   = creat_array(nx_all*nz_all, 0.0);
    Vz_2   = creat_array(nx_all*nz_all, 0.0);
    Vx_1   = creat_array(nx_all*nz_all, 0.0);
    Vx_2   = creat_array(nx_all*nz_all, 0.0);

    hTxx_1 = creat_array(nx_all*nz_all, 0.0);
    hTxx_2 = creat_array(nx_all*nz_all, 0.0);
    hTzz_1 = creat_array(nx_all*nz_all, 0.0);
    hTzz_2 = creat_array(nx_all*nz_all, 0.0);
    hTxz_1 = creat_array(nx_all*nz_all, 0.0);
    hTxz_2 = creat_array(nx_all*nz_all, 0.0);
    hVz_1  = creat_array(nx_all*nz_all, 0.0);
    hVz_2  = creat_array(nx_all*nz_all, 0.0);
    hVx_1  = creat_array(nx_all*nz_all, 0.0);
    hVx_2  = creat_array(nx_all*nz_all, 0.0);


/*======================================================
 * run
 *======================================================*/
    if (prepared_media == PREPARE_TTI) {

        ierr = write_model_TTI(save_model,
                c11_1, c13_1, c33_1,
                c15_1, c35_1, c55_1,
                c11_2, c13_2, c33_2,
                c15_2, c35_2, c55_2,
                B01, B10, nx_all, nz_all,
                nghost_x1, nghost_x2, nghost_z1, nghost_z2);

        ierr = elastic2d_lebedev(dx, dz, nx_all, nz_all, nt, dt,
                half_fd_stencil, spatial_difference_method, filter_method,
                xvec1, zvec1, xvec2, zvec2,
                c11_1, c13_1, c15_1, c33_1, c35_1, c55_1,
                c11_2, c13_2, c15_2, c33_2, c35_2, c55_2,
                B01, B10, src, source_impulse_method,
                boundary_type, boundary_layer_number,
                Txx_1 , Txx_2, Txz_1 , Txz_2, Tzz_1, Tzz_2,
                Vx_1  , Vx_2 , Vz_1, Vz_2,
                hTxx_1,hTxx_2, hTxz_1, hTxz_2, hTzz_1, hTzz_2,
                hVx_1 , hVx_2, hVz_1 , hVz_2,
                seismotype, nreceivers, xr, zr,
                NSTEP_BETWEEN_OUTPUT_SEISMOS,
                save_ASCII_seismograms, save_binary_seismograms,
                NSTEP_BETWEEN_OUTPUT_INFO,
                NSTEP_BETWEEN_OUTPUT_IMAGES,
                output_postscript_snapshot,
                imagetype_postscript, meshvect, modelvect, boundvect,
                sizemax_arrows, US_LETTER,
                output_wavefield_dumps, imagetype_wavefield_dumps,
                use_binary_for_wavefield_dumps);

    }
    else if (prepared_media == PREPARE_ISO) {
        ierr = elastic2d_staggered(dx, dz, nx_all, nz_all, nt, dt,
                    half_fd_stencil, spatial_difference_method,
                    xvec1, zvec1, xvec2, zvec2,
                    c11_1, c13_1, c11_1, c55_1,
                    B01, B10, src, source_impulse_method,
                    boundary_type, boundary_layer_number,
                    Txx_1, Txz_1, Tzz_1,
                    Vx_1,  Vz_1,
                    hTxx_1,hTxz_1, hTzz_1,
                    hVx_1, hVz_1,
                    seismotype, nreceivers, xr, zr,
                    NSTEP_BETWEEN_OUTPUT_SEISMOS,
                    save_ASCII_seismograms, save_binary_seismograms,
                    NSTEP_BETWEEN_OUTPUT_INFO,
                    NSTEP_BETWEEN_OUTPUT_IMAGES,
                    output_postscript_snapshot,
                    imagetype_postscript, meshvect, modelvect, boundvect,
                    sizemax_arrows, US_LETTER,
                    output_wavefield_dumps, imagetype_wavefield_dumps,
                    use_binary_for_wavefield_dumps);
    }
    else if (prepared_media == PREPARE_ORT) {
        ierr = elastic2d_staggered(dx, dz, nx_all, nz_all, nt, dt,
                    half_fd_stencil, spatial_difference_method,
                    xvec1, zvec1, xvec2, zvec2,
                    c11_1, c13_1, c33_1, c55_1,
                    B01, B10, src, source_impulse_method,
                    boundary_type, boundary_layer_number,
                    Txx_1, Txz_1, Tzz_1,
                    Vx_1,  Vz_1,
                    hTxx_1,hTxz_1, hTzz_1,
                    hVx_1, hVz_1,
                    seismotype, nreceivers, xr, zr,
                    NSTEP_BETWEEN_OUTPUT_SEISMOS,
                    save_ASCII_seismograms, save_binary_seismograms,
                    NSTEP_BETWEEN_OUTPUT_INFO,
                    NSTEP_BETWEEN_OUTPUT_IMAGES,
                    output_postscript_snapshot,
                    imagetype_postscript, meshvect, modelvect, boundvect,
                    sizemax_arrows, US_LETTER,
                    output_wavefield_dumps, imagetype_wavefield_dumps,
                    use_binary_for_wavefield_dumps);
    }

    free(xvec1); free(zvec1); free(xvec2); free(zvec2);
    free(c11_1); free(c13_1); free(c15_1); free(c33_1); free(c35_1);
    free(c11_2); free(c13_2); free(c15_2); free(c33_2); free(c35_2);
    free(lam2mu00); free(lam2mu11); free(lam11); free(lam00);
    free(mu00); free(mu11); free(rho_x); free(rho_z); free(rho01); free(rho10);
    free(B01  ); free(B10  ); free(c55_1); free(c55_2);
    free(Txx_1); free(Tzz_1); free(Txz_1); free(Vx_1); free(Vz_1);
    free(Txx_2); free(Tzz_2); free(Txz_2); free(Vx_2); free(Vz_2);
    free(hTxx_1); free(hTzz_1); free(hTxz_1); free(hVx_1); free(hVz_1);
    free(hTxx_2); free(hTzz_2); free(hTxz_2); free(hVx_2); free(hVz_2);
    free(xr); free(zr);
    free(src.stf_type_id); free(src.stf_timefactor); free(src.stf_freqfactor);
    free(src.xs); free(src.zs); free(src.Fx); free(src.Fz);
    free(src.Mxx); free(src.Mzz); free(src.Mxz);
    return 0;
}


float *creat_array(int n, float v)
{
    int i;
    float *array;
    array = (float*) malloc( n * sizeof(float));
    for (i = 0; i < n; i++ ) {
        array[i] = v;
    }
    return array;
}
