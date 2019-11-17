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
 *          11/2019: Add SSG code
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


    if (seismotype == SEIMO_VELOC) {
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
    else if (seismotype == SEIMO_STRESS) {
        write_stress_seismo(Txx_1, Txx_2,
                              Txz_1, Txz_2,
                              Tzz_1, Tzz_2,
                              nx,
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
    int *izvec1_rec, *izvec2_rec, *ixvec1_rec, *ixvec2_rec;
    float *DNx1, *DNx2, *DNz1, *DNz2;
    float Vx, Vz, damp;

    current_time = dt*(it+0.5);

    ixvec1_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    ixvec2_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    izvec1_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    izvec2_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));

    DNx1 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNz1 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNx2 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNz2 = (float*) malloc(2*NRECSAMPLE*sizeof(float));

    if(nreceiver <= 0)
        return;

    //-- for the time step
    if (it%NSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == nt-1) {
        //-- for every receiver
        for (ireceiver = 0; ireceiver < nreceiver; ireceiver++) {

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
                DNz2[i] = calculate_windowed_interp(izvec2_rec[i]-indx_rz_2, KAISERB, NRECSAMPLE);
            }


            Vx = 0.0;
            Vz = 0.0;

            for (i = 0; i < 2*NRECSAMPLE; i++){
                for (k = 0; k < 2*NRECSAMPLE; k++) {
                    // for Vx
                    indx = izvec1_rec[k] * nx + ixvec2_rec[i];
                    damp = DNx2[i] * DNz1[k]/4.0;
                    Vx += Vx_1[indx] * damp;
                    Vz += Vz_2[indx] * damp;

                    // for Vz
                    indx = izvec2_rec[k] * nx + ixvec1_rec[i];
                    damp = DNx1[i] * DNz2[k]/4.0;
                    Vz += Vz_1[indx] * damp;
                    Vx += Vx_2[indx] * damp;

                }
            }


            if (save_ASCII_seismograms) {

                //-- for Vx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vx.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6e\n", current_time, Vx);
                fclose(fp);

                //-- for Vz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vz.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6e\n", current_time, Vz);
                fclose(fp);
            }
            if (save_binary_seismograms) {
                //-- for Vx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vx.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Vx,sizeof(float),1,fp);
                fclose(fp);

                //-- for Vz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vz.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Vz,sizeof(float),1,fp);
                fclose(fp);
            }
        }

    }

    free(ixvec1_rec); free(izvec1_rec); free(ixvec2_rec); free(izvec2_rec);
    free(DNx1); free(DNx2); free(DNz1); free(DNz2);

    return 0;
}

int write_stress_seismo(float *Txx_1, float *Txx_2,
                          float *Txz_1, float *Txz_2,
                          float *Tzz_1, float *Tzz_2,
                          int nx,
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
    float Txx, Tzz, Txz, damp;

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
        if (it%NSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == nt-1) {

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

            //-- use average of two set of grids
            //
            Txx = 0.0;
            Txz = 0.0;
            Tzz = 0.0;

            for (i = 0; i < 2*NRECSAMPLE; i++){
                for (k = 0; k < 2*NRECSAMPLE; k++) {
                    // for Txx Tzz
                    indx = izvec1_rec[k] * nx + ixvec2_rec[i];
                    damp = DNx1[i] * DNz1[k]/4.0;
                    Txx += Txx_1[indx] * damp;
                    Tzz += Tzz_1[indx] * damp;
                    Txz += Txz_2[indx] * damp;

                    // for Txx Tzz
                    indx = izvec2_rec[k] * nx + ixvec1_rec[i];
                    damp = DNx2[i] * DNz2[k]/4.0;
                    Txx += Txx_2[indx] * damp;
                    Tzz += Tzz_2[indx] * damp;
                    Txz += Txz_1[indx] * damp;

                }
            }


            if (save_ASCII_seismograms) {
                //-- for Txx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txx.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6f\n", current_time, Txx);
                fclose(fp);

                //-- for Tzz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Tzz.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6f\n", current_time, Tzz );
                fclose(fp);

                //-- for Txz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txz.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6f\n", current_time, Txz );
                fclose(fp);
            }
            if (save_binary_seismograms) {
                //-- for Txx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txx.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Txx,sizeof(float),1,fp);
                fclose(fp);

                //-- for Tzz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Tzz.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Tzz,sizeof(float),1,fp);
                fclose(fp);

                //-- for Txz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txz.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Txz,sizeof(float),1,fp);
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
    if (abs(x) < 1e-5) {
        y = 1.0;
    }
    else if (1.0-x*x/r*r > 0.0) {
    /* Kai_sourceer finite impulse response filters. */
        y = first_modified_Bessel( 0, b*sqrt(1-x*x/r*r) ) / first_modified_Bessel(0, b);
    /* sinc */
        y = y * sin(PI*x) / (PI*x);
    }
    else {
        y = 0.0;
    }

    return y;
}


int write_seismograms_ssg(float *Vx, float *Vz,
                      float *Txx, float *Txz, float *Tzz,
                      float *xvec1, float *xvec2,
                      float *zvec1, float *zvec2,
                      int it, int nt, float dt,
                      int nx, int nz, float dx, float dz,
                      int nreceivers, float *xr, float *zr,
                      bool save_ASCII_seismograms, bool save_binary_seismograms,
                      int NSTEP_BETWEEN_OUTPUT_SEISMOS, int seismotype)
{


    if (seismotype == SEIMO_VELOC) {
        write_velocity_seismo_ssg(Vx, Vz, nx,
                              nreceivers, xr, zr,
                              it, nt, dt,
                              NSTEP_BETWEEN_OUTPUT_SEISMOS,
                              save_binary_seismograms,
                              save_ASCII_seismograms,
                              xvec1, xvec2,
                              zvec1, zvec2,
                              dx, dz);
    }
    else if (seismotype == SEIMO_STRESS) {
        write_stress_seismo_ssg(Txx, Txz, Tzz,
                                nx, nreceivers, xr, zr,
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

int write_velocity_seismo_ssg(float *Vx_1, float *Vz_1, int nx,
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
    int *izvec1_rec, *izvec2_rec, *ixvec1_rec, *ixvec2_rec;
    float *DNx1, *DNx2, *DNz1, *DNz2;
    float Vx, Vz, damp;

    current_time = dt*(it+0.5);

    ixvec1_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    ixvec2_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    izvec1_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));
    izvec2_rec = (int*)malloc(2*NRECSAMPLE*sizeof(int));

    DNx1 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNz1 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNx2 = (float*) malloc(2*NRECSAMPLE*sizeof(float));
    DNz2 = (float*) malloc(2*NRECSAMPLE*sizeof(float));

    if(nreceiver <= 0)
        return;

    //-- for the time step
    if (it%NSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == nt-1) {
        //-- for every receiver
        for (ireceiver = 0; ireceiver < nreceiver; ireceiver++) {

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
                DNz2[i] = calculate_windowed_interp(izvec2_rec[i]-indx_rz_2, KAISERB, NRECSAMPLE);
            }


            Vx = 0.0;
            Vz = 0.0;

            for (i = 0; i < 2*NRECSAMPLE; i++){
                for (k = 0; k < 2*NRECSAMPLE; k++) {
                    // for Vx
                    indx = izvec1_rec[k] * nx + ixvec2_rec[i];
                    damp = DNx2[i] * DNz1[k];
                    Vx += Vx_1[indx] * damp;

                    // for Vz
                    indx = izvec2_rec[k] * nx + ixvec1_rec[i];
                    damp = DNx1[i] * DNz2[k];
                    Vz += Vz_1[indx] * damp;

                }
            }


            if (save_ASCII_seismograms) {

                //-- for Vx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vx.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6e\n", current_time, Vx);
                fclose(fp);

                //-- for Vz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vz.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6e\n", current_time, Vz);
                fclose(fp);
            }
            if (save_binary_seismograms) {
                //-- for Vx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vx.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Vx,sizeof(float),1,fp);
                fclose(fp);

                //-- for Vz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Vz.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Vz,sizeof(float),1,fp);
                fclose(fp);
            }
        }

    }

    free(ixvec1_rec); free(izvec1_rec); free(ixvec2_rec); free(izvec2_rec);
    free(DNx1); free(DNx2); free(DNz1); free(DNz2);

    return 0;
}

int write_stress_seismo_ssg(float *Txx_1, float *Txz_1, float *Tzz_1,
                            int nx, int nreceiver, float *xr, float *zr,
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
    float Txx, Tzz, Txz, damp;

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
        if (it%NSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == nt-1) {

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

            //-- use average of two set of grids
            //
            Txx = 0.0;
            Txz = 0.0;
            Tzz = 0.0;

            for (i = 0; i < 2*NRECSAMPLE; i++){
                for (k = 0; k < 2*NRECSAMPLE; k++) {
                    // for Txx Tzz
                    indx = izvec1_rec[k] * nx + ixvec2_rec[i];
                    damp = DNx1[i] * DNz1[k];
                    Txx += Txx_1[indx] * damp;
                    Tzz += Tzz_1[indx] * damp;

                    // for Txx Tzz
                    indx = izvec2_rec[k] * nx + ixvec1_rec[i];
                    damp = DNx2[i] * DNz2[k];
                    Txz += Txz_1[indx] * damp;

                }
            }


            if (save_ASCII_seismograms) {
                //-- for Txx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txx.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6f\n", current_time, Txx);
                fclose(fp);

                //-- for Tzz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Tzz.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6f\n", current_time, Tzz );
                fclose(fp);

                //-- for Txz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txz.semd",ireceiver+1);
                fp = gfopen(file_name,"a");
                fprintf(fp, "%-20.6f%-20.6f\n", current_time, Txz );
                fclose(fp);
            }
            if (save_binary_seismograms) {
                //-- for Txx
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txx.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Txx,sizeof(float),1,fp);
                fclose(fp);

                //-- for Tzz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Tzz.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Tzz,sizeof(float),1,fp);
                fclose(fp);

                //-- for Txz
                snprintf(file_name, sizeof(file_name),"OUTPUT/S%04d.Txz.bin",it);
                fp = gfopen(file_name,"ab");
                fwrite(&current_time,sizeof(float),1,fp);
                fwrite(&Txz,sizeof(float),1,fp);
                fclose(fp);
            }
        }

    }

    free(ixvec1_rec); free(izvec1_rec); free(ixvec2_rec); free(izvec2_rec);
    free(DNx1); free(DNx2); free(DNz1); free(DNz2);

    return 0;
}
