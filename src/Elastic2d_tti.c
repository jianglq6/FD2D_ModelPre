
#include "elastic2d_src.h"
#include "Elastic2d_tti.h"
#include "elastic2d_lebedev.h"
#include "share_param.h"
#include "read_config_para.h"


int main()
{
    float vp = 4000.0, vs = 2700.0, rho = 2500.0;
    float lam2mu, mu, lam;
    lam2mu = vp*vp*rho;
    mu     = vs*vs*rho;
    lam    = lam2mu - 2*mu;
    fprintf(stdout, "lam2mu=%f, mu=%f, lam=%f\n",lam2mu, mu, lam);
    int ierr = 0, ib;
    char config_file[] = "./DATA/Par_file";
/*======================================================================
 *
 * read and initialize
 *
 * ====================================================================*/

/*============= effective media parameterization methond ==============*/
    int effective_para_method;
//--- filter
    int filter_method;
/*========================= time information ==========================*/
    int nt;
    float dt;
/*================== spatial difference information ===================*/
    int half_fd_stencil, spatial_difference_method;
/*========================= grid information ==========================*/
    int  nx, nz, ix, iz, indx;
    float  dx, dz, xmin, zmin;
/*========================= source information ========================*/
    struct Src src;
    int source_impulse_method, is;
/*======================== receiver info =============================*/
    int nreceivers, seismotype, NSTEP_BETWEEN_OUTPUT_SEISMOS;
    bool save_binary_seismograms, save_ASCII_seismograms;
    float *xr=NULL, *zr=NULL;
    //xr = (float*) malloc(1*sizeof(float));
    //zr = (float*) malloc(1*sizeof(float));
/*======================= boundary conditions ========================*/
    int boundary_type[4], boundary_layer_number[4];
/*=====================  velocity and density model ==================*/
    //-- size of model after add boundary layer
    int nx_all, nz_all;
    float xmin_all, zmin_all;
    float *xvec1, *xvec2, *zvec1, *zvec2;
    float *c11_1=NULL, *c13_1=NULL, *c15_1=NULL, *c33_1=NULL, *c35_1=NULL, *c55_1=NULL;
    float *c11_2=NULL, *c13_2=NULL, *c15_2=NULL, *c33_2=NULL, *c35_2=NULL, *c55_2=NULL;
    float *B01=NULL, *B10=NULL;
/*=================== display parameters ============================*/
    int NSTEP_BETWEEN_OUTPUT_INFO;
/*============== modvies snapshot information =======================*/
    int NSTEP_BETWEEN_OUTPUT_IMAGES;
    bool output_wavefield_dumps, use_binary_for_wavefield_dumps;
    int imagetype_wavefield_dumps=1, imagetype_postscript=2;
    bool meshvect, modelvect, boundvect, US_LETTER, output_postscript_snapshot;
    float sizemax_arrows=1.0;
/*======================= alloc for Lebedev =========================*/
    float *Txx_1=NULL, *Txx_2=NULL, *Tzz_1=NULL, *Tzz_2=NULL, *Txz_1=NULL, *Txz_2=NULL;
    float *Vx_1=NULL,  *Vx_2=NULL,  *Vz_1=NULL,  *Vz_2=NULL;
    float *hTxx_1=NULL, *hTxx_2=NULL, *hTzz_1=NULL, *hTzz_2=NULL, *hTxz_1=NULL, *hTxz_2=NULL;
    float *hVx_1=NULL,  *hVx_2=NULL,  *hVz_1=NULL,  *hVz_2=NULL;

    /*==== get the info and print ====*/
    ierr = get_config_info(config_file, &nt, &dt,
                           &half_fd_stencil, &spatial_difference_method,
                           &xmin, &dx, &nx, &zmin, &dz, &nz, &filter_method,
                           &source_impulse_method, &src,
                           &seismotype, &NSTEP_BETWEEN_OUTPUT_SEISMOS,
                           &save_ASCII_seismograms, &save_binary_seismograms,
                           &nreceivers, &xr, &zr,
                            boundary_type, boundary_layer_number,
                           &NSTEP_BETWEEN_OUTPUT_INFO,
                           &NSTEP_BETWEEN_OUTPUT_IMAGES, &output_postscript_snapshot,
                           &imagetype_postscript, &meshvect, &modelvect, &boundvect,
                           &sizemax_arrows, &US_LETTER,
                           &output_wavefield_dumps, &imagetype_wavefield_dumps,
                           &use_binary_for_wavefield_dumps);

    /* prepare for boundary condition */
    for (ib = 0; ib < 4; ib++){
        if (boundary_type[ib] == 0)
            boundary_layer_number[ib] = 0;
    }

    /* reconstruct the coordinate */
    xmin_all = xmin - (boundary_layer_number[0]+half_fd_stencil) * dx;
    zmin_all = zmin - (boundary_layer_number[2]+half_fd_stencil) * dz;
    nx_all = nx + boundary_layer_number[0] + boundary_layer_number[1] + half_fd_stencil*2;
    nz_all = nz + boundary_layer_number[2] + boundary_layer_number[2] + half_fd_stencil*2;

    zvec1 = (float*) malloc(nz_all*sizeof(float));
    zvec2 = (float*) malloc(nz_all*sizeof(float));
    xvec1 = (float*) malloc(nx_all*sizeof(float));
    xvec2 = (float*) malloc(nx_all*sizeof(float));

    /* alloc the material array */
    c11_1 = creat_component_array(nx_all*nz_all);
    c13_1 = creat_component_array(nx_all*nz_all);
    c15_1 = creat_component_array(nx_all*nz_all);
    c33_1 = creat_component_array(nx_all*nz_all);
    c35_1 = creat_component_array(nx_all*nz_all);
    c55_1 = creat_component_array(nx_all*nz_all);

    c11_2 = creat_component_array(nx_all*nz_all);
    c13_2 = creat_component_array(nx_all*nz_all);
    c15_2 = creat_component_array(nx_all*nz_all);
    c33_2 = creat_component_array(nx_all*nz_all);
    c35_2 = creat_component_array(nx_all*nz_all);
    c55_2 = creat_component_array(nx_all*nz_all);

    B10   = creat_component_array(nx_all*nz_all);
    B01   = creat_component_array(nx_all*nz_all);

    for (ix=0; ix<nx_all; ix++) {
        xvec1[ix] = xmin_all      + dx*ix;
        xvec2[ix] = xmin_all+dx/2 + dx*ix;
    }

    for (iz=0; iz<nz_all; iz++) {
        zvec1[iz] = zmin_all      + dz*iz;
        zvec2[iz] = zmin_all+dz/2 + dz*iz;
    }


/*!!! just for test */
    for (ix = 0; ix<nx_all; ix++) {
        for (iz=0; iz<nz_all; iz++) {
            indx = iz*nx_all + ix;
        //-- assign
            c11_1[indx] = lam2mu;
            c11_2[indx] = lam2mu;
            c13_1[indx] = lam;
            c13_2[indx] = lam;
            c15_1[indx] = 0.0;
            c15_2[indx] = 0.0;
            c33_1[indx] = lam2mu;
            c33_2[indx] = lam2mu;
            c35_1[indx] = 0.0;
            c35_2[indx] = 0.0;
            c55_1[indx] = mu;
            c55_2[indx] = mu;
            B01[indx]   = 1/rho;
            B10[indx]   = 1/rho;

        }
    }



    /* alloc the components */
    Txx_1  = creat_component_array(nx_all*nz_all);
    Txx_2  = creat_component_array(nx_all*nz_all);
    Tzz_1  = creat_component_array(nx_all*nz_all);
    Tzz_2  = creat_component_array(nx_all*nz_all);
    Txz_1  = creat_component_array(nx_all*nz_all);
    Txz_2  = creat_component_array(nx_all*nz_all);
    Vz_1   = creat_component_array(nx_all*nz_all);
    Vz_2   = creat_component_array(nx_all*nz_all);
    Vx_1   = creat_component_array(nx_all*nz_all);
    Vx_2   = creat_component_array(nx_all*nz_all);

    hTxx_1 = creat_component_array(nx_all*nz_all);
    hTxx_2 = creat_component_array(nx_all*nz_all);
    hTzz_1 = creat_component_array(nx_all*nz_all);
    hTzz_2 = creat_component_array(nx_all*nz_all);
    hTxz_1 = creat_component_array(nx_all*nz_all);
    hTxz_2 = creat_component_array(nx_all*nz_all);
    hVz_1  = creat_component_array(nx_all*nz_all);
    hVz_2  = creat_component_array(nx_all*nz_all);
    hVx_1  = creat_component_array(nx_all*nz_all);
    hVx_2  = creat_component_array(nx_all*nz_all);

/*======================================================
 * run
 *======================================================*/

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


    free(xvec1); free(zvec1); free(xvec2); free(zvec2);
    free(c11_1); free(c13_1); free(c15_1); free(c33_1); free(c35_1);
    free(c11_2); free(c13_2); free(c15_2); free(c33_2); free(c35_2);
    free(B01  ); free(B10  );
    free(Txx_1); free(Tzz_1); free(Txz_1); free(Vx_1); free(Vz_1);
    free(Txx_2); free(Tzz_2); free(Txz_2); free(Vx_2); free(Vz_2);
    free(hTxx_1); free(hTzz_1); free(hTxz_1); free(hVx_1); free(hVz_1);
    free(hTxx_2); free(hTzz_2); free(hTxz_2); free(hVx_2); free(hVz_2);
    free(xr); free(zr);
    //free(src.stf_type_id); free(src.stf_timefactor); free(src.stf_freqfactor);
    //free(src.xs); free(src.zs); free(src.Fx); free(src.Fz);
    //free(src.Mxx); free(src.Mzz); free(src.Mxz);
    return 0;
}

float *creat_component_array(int n)
{
    int i;
    float *array;
    array = (float*) malloc( n * sizeof(float));
    for (i = 0; i < n; i++ ) {
        array[i] = 0.0;
    }
    return array;
}
