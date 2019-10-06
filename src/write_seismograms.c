/***************************************************************************
 *
 * This function is used to write seismograms to text or binary files
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 09/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/

#include "read_config_para.h"
#include "elastic2d_math.h"
#include "write_seismograms.h"

int write_seismograms(float *Vx_1, float *Vx_2,
                      float *Vz_1, float *Vz_2,
                      float *Txx_1, float *Tzz_1,
                      float *Txz_1, float *Txx_2,
                      float *Tzz_2, float *Txz_2,
                      float *xvec1, float *xvec2,
                      float *zvec1, float *zvec2,
                      int it, int nt, float dt,
                      int nx, int nz, float dx, float dz,
                      int nreceivers, float *xr, float *zr,
                      bool save_ASCII_seismograms, bool save_binary_seismograms,
                      int NSTEP_BETWEEN_OUTPUT_SEISMOS, int seismotype)
{

    if (seismotype == 1) {
        write_velocity_seismo(Vx_1, Vx_2,
                              Vz_1, Vz_2, nx,
                              nreceivers, xr, zr,
                              it, nt, dt,
                              NSTEP_BETWEEN_OUTPUT_SEISMOS,
                              save_binary_seismograms,
                              save_ASCII_seismograms,
                              xvec1, xvec2,
                              zvec1, zvec2,
                              dx, dz);
    }


    return 0;

}

int write_velocity_seismo(float *Vx_1, float *Vx_2,
                          float *Vz_1, float *Vz_2, int nx,
                          int nreceiver, float *xr, float *zr,
                          int it, int nt, float dt,
                          int NSTEP_BETWEEN_OUTPUT_SEISMOS,
                          bool save_binary_seismograms,
                          bool save_ASCII_seismograms,
                          float *xvec1, float *xvec2,
                          float *zvec1, float *zvec2,
                          float dx, float dz)
{
    int ireceiver, i, k, indx;
    FILE *fp;
    char file_name[MAX_BUF_LEN];
    float current_time, indx_rx_1, indx_rz_1, indx_rx_2, indx_rz_2;
    int *izvec1_rec=NULL, *izvec2_rec=NULL, *ixvec1_rec=NULL, *ixvec2_rec=NULL;
    float *DNx1=NULL, *DNx2=NULL, *DNz1=NULL, *DNz2=NULL;
    float Vx, Vz, damp;

    current_time = dt*it;

    ixvec1_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    ixvec2_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    izvec1_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    izvec2_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));

    DNx1 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNz1 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNx2 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNz2 = (float*) malloc(2*NRECSAMPLE*sizeof(float));

    //-- for every receiver
    for (ireceiver = 0; ireceiver < nreceiver; ireceiver++) {
        //-- for the time step
        if (it%NSTEP_BETWEEN_OUTPUT_SEISMOS == 0 && it == nt-1) {

            indx_rx_1 = (xr[ireceiver] - xvec1[0])/dx;
            indx_rz_1 = (zr[ireceiver] - zvec1[0])/dz;
            indx_rx_2 = (xr[ireceiver] - xvec2[0])/dx;
            indx_rz_2 = (zr[ireceiver] - zvec2[0])/dz;

            /* calculate index of receiver distribution */
            for (i = 0; i < NRECSAMPLE; i++) {
                ixvec1_rec[NRECSAMPLE-i-1] = floor(indx_rx_1) - i;
                ixvec1_rec[NRECSAMPLE+i  ] =  ceil(indx_rx_1) + i;
                izvec1_rec[NRECSAMPLE-i-1] = floor(indx_rz_1) - i;
                izvec1_rec[NRECSAMPLE+i  ] =  ceil(indx_rz_1) + i;

                ixvec2_rec[NRECSAMPLE-i-1] = floor(indx_rx_2) - i;
                ixvec2_rec[NRECSAMPLE+i  ] =  ceil(indx_rx_2) + i;
                izvec2_rec[NRECSAMPLE-i-1] = floor(indx_rz_2) - i;
                izvec2_rec[NRECSAMPLE+i  ] =  ceil(indx_rz_2) + i;
            }

            /* calculate seismo with distribution receiver */
            for (i = 0; i < 2*NRECSAMPLE; i++) {
                DNx1[i] = calculate_windowed_interp(ixvec1_rec[i]-indx_rx_1, KAISERB, NRECSAMPLE);
                DNx2[i] = calculate_windowed_interp(ixvec2_rec[i]-indx_rx_2, KAISERB, NRECSAMPLE);
                DNz1[i] = calculate_windowed_interp(izvec1_rec[i]-indx_rz_1, KAISERB, NRECSAMPLE);
                DNz2[i] = calculate_windowed_interp(izvec1_rec[i]-indx_rz_2, KAISERB, NRECSAMPLE);
            }

        //-- TODO! plan 1: just use the  first set of the grid
            Vx = 0.0;
            Vz = 0.0;

            for (i = 0; i < 2*NRECSAMPLE; i++){
                for (k = 0; k < 2*NRECSAMPLE; k++) {
                    // for Vx
                    indx = izvec1_rec[k] * nx + ixvec2_rec[i];
                    damp = DNx2[i] * DNz1[k];
                    Vx = Vx + Vx_1[indx] * damp;

                    // for Vz
                    indx = izvec2_rec[k] * nx + ixvec1_rec[i];
                    damp = DNx1[i] * DNz2[k];
                    Vz = Vz + Vz_1[indx] * damp;

                }
            }


            if (save_ASCII_seismograms) {
                //-- for Vx
                snprintf(file_name, sizeof(file_name),"S%04d.Vx.semd",ireceiver+1);
                fp = fopen(file_name,"a");
                fprintf(fp, "%4.6f    %1.12e\n", current_time, Vx);
                fclose(fp);

                //-- for Vx
                snprintf(file_name, sizeof(file_name),"S%04d.Vz.semd",ireceiver+1);
                fp = fopen(file_name,"a");
                fprintf(fp, "%4.6f    %1.12e\n", current_time, Vz );
                fclose(fp);
            }
        }

    }

    free(ixvec1_rec); free(izvec1_rec); free(ixvec2_rec); free(izvec2_rec);
    free(DNx1); free(DNx2); free(DNz1); free(DNz2);

    return 0;
}


/* Dicrete receiver function                                               *
 *  ref: Hicks, 2002, Arbitrary source and receiver positioning in       *
 *        finite-difference schemes using Kaier windowwd sinc function  */
float calculate_windowed_interp(float x, float b, float r)
{
    float y = 0.0;
    if (1.0-x*x/r*r > 0.0) {
    /* interpolation coefficient */
        y = first_modified_Bessel( 0, b*sqrt(1-x*x/r*r) ) / first_modified_Bessel(0, b);
    /* sinc */
        y = y * sin(PI*x) / (PI*x);
    }
    else
        y = 0.0;

    return y;
}


