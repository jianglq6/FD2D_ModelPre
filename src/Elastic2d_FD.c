
#include "utility_vector.h"
#include "elastic2d_src.h"
#include "elastic2d_lebedev.h"
#include "elastic2d_staggered.h"
#include "share_param.h"
#include "read_config_para.h"
#include "pre_model_prepare.h"

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
    int boundary_type, boundary_layer_number[4], ib;

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
    float *Txx00=NULL, *Txx11=NULL, *Tzz00=NULL, *Tzz11=NULL, *Txz11=NULL, *Txz00=NULL;
    float *Vx01=NULL,  *Vx10=NULL,  *Vz10=NULL,  *Vz01=NULL;
    float *hTxx00=NULL, *hTxx11=NULL, *hTzz00=NULL, *hTzz11=NULL, *hTxz11=NULL, *hTxz00=NULL;
    float *hVx01=NULL,  *hVx10=NULL,  *hVz10=NULL,  *hVz01=NULL;
    float *DxTxx01=NULL, *DzTxz01=NULL, *DxTxz01=NULL, *DzTzz01=NULL;
    float *DxTxx10=NULL, *DzTxz10=NULL, *DxTxz10=NULL, *DzTzz10=NULL;
    float *DxVx00=NULL, *DzVz00=NULL, *DzVx00=NULL, *DxVz00=NULL; /* partial difference centered of (i,j) */
    float *DxVx11=NULL, *DzVz11=NULL, *DzVx11=NULL, *DxVz11=NULL; /* partial difference centered of (i+1/2,j+1/2) */

/*===================== get the info and print =======================*/
    ierr = get_config_info(config_file, &nt, &dt,
                           &half_fd_stencil, &spatial_difference_method,
                           &xmin, &dx, &nx, &zmin, &dz, &nz, &filter_method,
                           &source_impulse_method, &src,
                           &seismotype, &NSTEP_BETWEEN_OUTPUT_SEISMOS,
                           &save_ASCII_seismograms, &save_binary_seismograms,
                           &nreceivers, &xr, &zr,
                           &boundary_type, boundary_layer_number,
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
    c11_1 = creat_float_array(nx_all*nz_all, -1.0);
    c13_1 = creat_float_array(nx_all*nz_all, -1.0);
    c15_1 = creat_float_array(nx_all*nz_all,  0.0);
    c33_1 = creat_float_array(nx_all*nz_all, -1.0);
    c35_1 = creat_float_array(nx_all*nz_all,  0.0);
    c55_1 = creat_float_array(nx_all*nz_all, -1.0);

    c11_2 = creat_float_array(nx_all*nz_all, -1.0);
    c13_2 = creat_float_array(nx_all*nz_all, -1.0);
    c15_2 = creat_float_array(nx_all*nz_all,  0.0);
    c33_2 = creat_float_array(nx_all*nz_all, -1.0);
    c35_2 = creat_float_array(nx_all*nz_all,  0.0);
    c55_2 = creat_float_array(nx_all*nz_all, -1.0);

    B10   = creat_float_array(nx_all*nz_all, 1.0);
    B01   = creat_float_array(nx_all*nz_all, 1.0);

    lam2mu00 = creat_float_array(nx_all*nz_all, -1.0);
    lam2mu11 = creat_float_array(nx_all*nz_all, -1.0);
    lam11    = creat_float_array(nx_all*nz_all, -1.0);
    lam00    = creat_float_array(nx_all*nz_all, -1.0);
    mu00     = creat_float_array(nx_all*nz_all, -1.0);
    mu11     = creat_float_array(nx_all*nz_all, -1.0);
    rho_x    = creat_float_array(nx_all*nz_all, -1.0);
    rho_z    = creat_float_array(nx_all*nz_all, -1.0);
    rho01    = creat_float_array(nx_all*nz_all, -1.0);
    rho10    = creat_float_array(nx_all*nz_all, -1.0);

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
                xvec1, zvec1, xvec2, zvec2, nx_all, nz_all, dx, dz,
                c11_1, c13_1, c15_1, c33_1, c35_1, c55_1,
                c11_2, c13_2, c15_2, c33_2, c35_2, c55_2,
                B01, B10,
                lam2mu00, lam2mu11, lam00, lam11, mu00, mu11,
                rho01, rho10, rho_x, rho_z);
    }



    /* alloc the components */
    Txx00  = creat_float_array(nx_all*nz_all, 0.0);
    Txx11  = creat_float_array(nx_all*nz_all, 0.0);
    Tzz00  = creat_float_array(nx_all*nz_all, 0.0);
    Tzz11  = creat_float_array(nx_all*nz_all, 0.0);
    Txz11  = creat_float_array(nx_all*nz_all, 0.0);
    Txz00  = creat_float_array(nx_all*nz_all, 0.0);
    Vz10   = creat_float_array(nx_all*nz_all, 0.0);
    Vz01   = creat_float_array(nx_all*nz_all, 0.0);
    Vx01   = creat_float_array(nx_all*nz_all, 0.0);
    Vx10   = creat_float_array(nx_all*nz_all, 0.0);

    hTxx00 = creat_float_array(nx_all*nz_all, 0.0);
    hTxx11 = creat_float_array(nx_all*nz_all, 0.0);
    hTzz00 = creat_float_array(nx_all*nz_all, 0.0);
    hTzz11 = creat_float_array(nx_all*nz_all, 0.0);
    hTxz11 = creat_float_array(nx_all*nz_all, 0.0);
    hTxz00 = creat_float_array(nx_all*nz_all, 0.0);
    hVz10  = creat_float_array(nx_all*nz_all, 0.0);
    hVz01  = creat_float_array(nx_all*nz_all, 0.0);
    hVx01  = creat_float_array(nx_all*nz_all, 0.0);
    hVx10  = creat_float_array(nx_all*nz_all, 0.0);

    DxVz00 = creat_float_array(nx_all*nz_all, 0.0);
    DxVz11 = creat_float_array(nx_all*nz_all, 0.0);
    DzVx00 = creat_float_array(nx_all*nz_all, 0.0);
    DzVx11 = creat_float_array(nx_all*nz_all, 0.0);
    DzVz00 = creat_float_array(nx_all*nz_all, 0.0);
    DzVz11 = creat_float_array(nx_all*nz_all, 0.0);
    DxVx00 = creat_float_array(nx_all*nz_all, 0.0);
    DxVx11 = creat_float_array(nx_all*nz_all, 0.0);

    DxTxx01 = creat_float_array(nx_all*nz_all, 0.0);
    DxTxx10 = creat_float_array(nx_all*nz_all, 0.0);
    DzTxz01 = creat_float_array(nx_all*nz_all, 0.0);
    DzTxz10 = creat_float_array(nx_all*nz_all, 0.0);
    DxTxz01 = creat_float_array(nx_all*nz_all, 0.0);
    DxTxz10 = creat_float_array(nx_all*nz_all, 0.0);
    DzTzz01 = creat_float_array(nx_all*nz_all, 0.0);
    DzTzz10 = creat_float_array(nx_all*nz_all, 0.0);



/*======================================================
 * run
 *======================================================*/
    if (prepared_media == PREPARE_TTI) {

        ierr = write_model_TTI(save_model,
                c11_1, c13_1, c33_1,
                c15_1, c35_1, c55_1,
                c11_2, c13_2, c33_2,
                c15_2, c35_2, c55_2,
                B01, B10, nx_all, nz_all);

        ierr = elastic2d_lebedev(dx, dz, nx_all, nz_all, nt, dt,
                half_fd_stencil, spatial_difference_method, filter_method,
                xvec1, zvec1, xvec2, zvec2,
                c11_1, c13_1, c15_1, c33_1, c35_1, c55_1,
                c11_2, c13_2, c15_2, c33_2, c35_2, c55_2,
                B01, B10, src, source_impulse_method,
                boundary_type, boundary_layer_number,
                Txx00 , Txx11, Txz11 , Txz00, Tzz00, Tzz11,
                Vx01  , Vx10 , Vz10, Vz01,
                hTxx00,hTxx11, hTxz11, hTxz00, hTzz00, hTzz11,
                hVx01 , hVx10, hVz10 , hVz01,
                DxVx00, DzVz00, DzVx00, DxVz00,
                DxVx11, DzVz11, DzVx11, DxVz11,
                DxTxx01, DzTxz01, DxTxz01, DzTzz01,
                DxTxx10, DzTxz10, DxTxz10, DzTzz10,
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
                    Txx00, Txz11, Tzz00,
                    Vx01,  Vz10,
                    hTxx00,hTxz11, hTzz00,
                    hVx01, hVz10,
                    DxVx00, DzVz00, DzVx11, DxVz11,
                    DxTxx01,DzTxz01,DxTxz10,DzTzz10,
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
                    Txx00, Txz11, Tzz00,
                    Vx01,  Vz10,
                    hTxx00,hTxz11, hTzz00,
                    hVx01, hVz10,
                    DxVx00, DzVz00, DzVx11, DxVz11,
                    DxTxx01,DzTxz01,DxTxz10,DzTzz10,
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
    free(Txx00); free(Tzz00); free(Txz11); free(Vx01); free(Vz10);
    free(Txx11); free(Tzz11); free(Txz00); free(Vx10); free(Vz01);
    free(hTxx00); free(hTzz00); free(hTxz11); free(hVx01); free(hVz10);
    free(hTxx11); free(hTzz11); free(hTxz00); free(hVx10); free(hVz01);
    free(DxTxx01); free(DzTxz01); free(DxTxz01); free(DzTzz01);
    free(DxTxx10); free(DzTxz10); free(DxTxz10); free(DzTzz10);
    free(DxVx00); free(DzVz00); free(DzVx00); free(DxVz00);
    free(DxVx11); free(DzVz11); free(DzVx11); free(DxVz11);
    free(xr); free(zr);
    free(src.stf_type_id); free(src.stf_timefactor); free(src.stf_freqfactor);
    free(src.xs); free(src.zs); free(src.Fx); free(src.Fz);
    free(src.Mxx); free(src.Mzz); free(src.Mxz);
    return 0;
}


