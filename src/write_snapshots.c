/***************************************************************************
 *
 * This function is used to write snapshots to ./OUTPUT file
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 10/2019: Original version created by Luqian Jiang
 *
 ***************************************************************************/
#include "write_snapshots.h"
#include "read_config_para.h"

int write_snapshots(float *Vx_1,  float *Vx_2,
                    float *Vz_1,  float *Vz_2,
                    float *Txx_1, float *Tzz_1,
                    float *Txz_1, float *Txx_2,
                    float *Tzz_2, float *Txz_2,
                    float *xvec1, float *xvec2,
                    float *zvec1, float *zvec2,
                    int it, int nt, float dt,   // it = it+1
                    int nx, int nz, float dx, float dz, int half_fd_stencil,
                    int *boundary_layer_number,
                    int NSTEP_BETWEEN_OUTPUT_IMAGES,
                    bool output_postscript_snapshot,
                    int imagetype_postscript, bool meshvect, bool modelvect, bool boundvect,
                    float sizemax_arrows, bool US_LETTER,
                    bool output_wavefield_dumps, int imagetype_wavefield_dumps,
                    bool use_binary_for_wavefield_dumps)
{

    if (it%NSTEP_BETWEEN_OUTPUT_IMAGES == 0 || it==5) {

        //if(output_postscript_snapshot) {
        //    write_postscript(Vx_1 , Vx_2,
        //                     Vz_1 , Vz_2,
        //                     Txx_1, Tzz_1,
        //                     Txz_1, Txx_2,
        //                     Tzz_2, Txz_2,
        //                     xvec1, xvec2,
        //                     zvec1, zvec2,
        //                     it, dt,
        //                     nx, nz, dx, dz,
        //                     boundary_layer_number,
        //                     NSTEP_BETWEEN_OUTPUT_IMAGES,
        //                     imagetype_postscript, meshvect, modelvect, boundvect,
        //                     sizemax_arrows, US_LETTER);
        //}

        if (output_wavefield_dumps) {
            write_wavefield_dumps(Vx_1 , Vx_2,
                                  Vz_1 , Vz_2,
                                  Txx_1, Tzz_1,
                                  Txz_1, Txx_2,
                                  Tzz_2, Txz_2,
                                  xvec1, xvec2,
                                  zvec1, zvec2,
                                  it, nt, dt,
                                  nx, nz, dx, dz,
                                  boundary_layer_number,
                                  imagetype_wavefield_dumps,
                                  use_binary_for_wavefield_dumps);
        }

    }

    return 0;

}

int write_wavefield_dumps(float *Vx_1,  float *Vx_2,
                          float *Vz_1,  float *Vz_2,
                          float *Txx_1, float *Tzz_1,
                          float *Txz_1, float *Txx_2,
                          float *Tzz_2, float *Txz_2,
                          float *xvec1, float *xvec2,
                          float *zvec1, float *zvec2,
                          int it, int nt, float dt,   // it = it+1
                          int nx, int nz, float dx, float dz,
                          int *boundary_layer_number,
                          int imagetype_wavefield_dumps,
                          bool use_binary_for_wavefield_dumps)
{
    FILE *fp_xvec1=NULL, *fp_xvec2=NULL, *fp_zvec1=NULL, *fp_zvec2=NULL;
    FILE *fp_field=NULL;
    int ix, iz;
    char fn[1024];

    if(it==5){
        if (use_binary_for_wavefield_dumps){
            fp_xvec1 = gfopen("OUTPUT/wavefield_xvec1_for_dumps.bin","wb");
            fp_xvec2 = gfopen("OUTPUT/wavefield_xvec2_for_dumps.bin","wb");
            fp_zvec1 = gfopen("OUTPUT/wavefield_zvec1_for_dumps.bin","wb");
            fp_zvec2 = gfopen("OUTPUT/wavefield_zvec2_for_dumps.bin","wb");
            for (ix = 0; ix < nx; ix++){
                fwrite(&xvec1[ix], sizeof(float), 1, fp_xvec1);
                fwrite(&xvec2[ix], sizeof(float), 1, fp_xvec2);
            }
            for (iz = 0; iz < nz; iz++){
                fwrite(&zvec1[iz], sizeof(float), 1, fp_zvec1);
                fwrite(&zvec2[iz], sizeof(float), 1, fp_zvec2);
            }
        }
        else {
            fp_xvec1 = gfopen("OUTPUT/wavefield_xvec1_for_dumps.txt","w");
            fp_xvec2 = gfopen("OUTPUT/wavefield_xvec2_for_dumps.txt","w");
            fp_zvec1 = gfopen("OUTPUT/wavefield_zvec1_for_dumps.txt","w");
            fp_zvec2 = gfopen("OUTPUT/wavefield_zvec2_for_dumps.txt","w");
            for (ix = 0; ix < nx; ix++){
                fprintf(fp_xvec1, "%20.6f\n", xvec1[ix]);
                fprintf(fp_xvec2, "%20.6f\n", xvec2[ix]);
            }
            for (iz = 0; iz < nz; iz++){
                fprintf(fp_zvec1, "%20.6f\n", zvec1[iz]);
                fprintf(fp_zvec2, "%20.6f\n", zvec2[iz]);
            }
        }


        fclose(fp_xvec1); fclose(fp_xvec2);
        fclose(fp_zvec1); fclose(fp_zvec2);

    }

    /* write the velocity components */
    if (imagetype_wavefield_dumps==SNAP_VELOC) {
        if (use_binary_for_wavefield_dumps) {
            /*For Vx1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vx1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Vx_1[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
            /* For vx2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vx2.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Vx_2[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

            /*For Vz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vz1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Vz_1[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
            /* For vx2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vz2.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Vz_2[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
        }
        /* write ASCII version */
        else {
            /*For Vx1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vx1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Vx_1[iz*nx+ix]);
                fprintf(fp_field,"\n");
            }
            fclose(fp_field);
            /* For Vx2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vx2.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Vx_2[iz*nx+ix]);
                fprintf(fp_field,"\n");
            }
            fclose(fp_field);

            /*For Vz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vz1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Vz_1[iz*nx+ix]);
                fprintf(fp_field,"\n");
            }
            fclose(fp_field);
            /* For Vz2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vz2.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Vz_2[iz*nx+ix]);
                fprintf(fp_field,"\n");
            }
            fclose(fp_field);

        }
    }

 /* write the velocity components */
    if (imagetype_wavefield_dumps==SNAP_STRESS) {
        if (use_binary_for_wavefield_dumps) {
            /*For Txx1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txx1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Txx_1[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
            /* For Txx2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txx2.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Txx_2[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

            /*For Tzz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Tzz1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Tzz_1[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
            /* For Tzz2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Tzz2.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Tzz_2[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

            /*For Txz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txz1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Txz_1[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
            /* For Txz2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txz2.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Txz_2[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);
        }
        /* write ASCII version */
        else {
            /*For Txx1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txx1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Txx_1[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);
            /* For Txx2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txx2.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Txx_2[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);

            /*For Tzz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Tzz1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Tzz_1[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);
            /* For Tzz2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Tzz2.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Tzz_2[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);

            /*For Txz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txz1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Txz_1[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);
            /* For Txz2 */
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txz2.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Txz_2[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);

        }
    }


    return 0;
}

int write_snapshots_ssg(float *Vx, float *Vz,
                        float *Txx, float *Tzz, float *Txz,
                        float *xvec1, float *xvec2,
                        float *zvec1, float *zvec2,
                        int it, int nt, float dt,   // it = it+1
                        int nx, int nz, float dx, float dz, int half_fd_stencil,
                        int *boundary_layer_number,
                        int NSTEP_BETWEEN_OUTPUT_IMAGES,
                        bool output_postscript_snapshot,
                        int imagetype_postscript, bool meshvect, bool modelvect, bool boundvect,
                        float sizemax_arrows, bool US_LETTER,
                        bool output_wavefield_dumps, int imagetype_wavefield_dumps,
                        bool use_binary_for_wavefield_dumps)
{

    if (it%NSTEP_BETWEEN_OUTPUT_IMAGES == 0 || it==5) {

        //if(output_postscript_snapshot) {
        //    write_postscript(Vx_1 , Vx_2,
        //                     Vz_1 , Vz_2,
        //                     Txx_1, Tzz_1,
        //                     Txz_1, Txx_2,
        //                     Tzz_2, Txz_2,
        //                     xvec1, xvec2,
        //                     zvec1, zvec2,
        //                     it, dt,
        //                     nx, nz, dx, dz,
        //                     boundary_layer_number,
        //                     NSTEP_BETWEEN_OUTPUT_IMAGES,
        //                     imagetype_postscript, meshvect, modelvect, boundvect,
        //                     sizemax_arrows, US_LETTER);
        //}

        if (output_wavefield_dumps) {
            write_wavefield_dumps_ssg(Vx , Vz,
                                  Txx, Tzz, Txz,
                                  xvec1, xvec2,
                                  zvec1, zvec2,
                                  it, nt, dt,
                                  nx, nz, dx, dz,
                                  boundary_layer_number,
                                  imagetype_wavefield_dumps,
                                  use_binary_for_wavefield_dumps);
        }

    }

    return 0;

}

int write_wavefield_dumps_ssg(float *Vx, float *Vz,
                              float *Txx, float *Tzz, float *Txz,
                              float *xvec1, float *xvec2,
                              float *zvec1, float *zvec2,
                              int it, int nt, float dt,   // it = it+1
                              int nx, int nz, float dx, float dz,
                              int *boundary_layer_number,
                              int imagetype_wavefield_dumps,
                              bool use_binary_for_wavefield_dumps)
{
    FILE *fp_xvec1=NULL, *fp_xvec2=NULL, *fp_zvec1=NULL, *fp_zvec2=NULL;
    FILE *fp_field=NULL;
    int ix, iz;
    char fn[1024];

    if(it==5){
        if (use_binary_for_wavefield_dumps){
            fp_xvec1 = gfopen("OUTPUT/wavefield_xvec1_for_dumps.bin","wb");
            fp_xvec2 = gfopen("OUTPUT/wavefield_xvec2_for_dumps.bin","wb");
            fp_zvec1 = gfopen("OUTPUT/wavefield_zvec1_for_dumps.bin","wb");
            fp_zvec2 = gfopen("OUTPUT/wavefield_zvec2_for_dumps.bin","wb");
            for (ix = 0; ix < nx; ix++){
                fwrite(&xvec1[ix], sizeof(float), 1, fp_xvec1);
                fwrite(&xvec2[ix], sizeof(float), 1, fp_xvec2);
            }
            for (iz = 0; iz < nz; iz++){
                fwrite(&zvec1[iz], sizeof(float), 1, fp_zvec1);
                fwrite(&zvec2[iz], sizeof(float), 1, fp_zvec2);
            }
        }
        else {
            fp_xvec1 = gfopen("OUTPUT/wavefield_xvec1_for_dumps.txt","w");
            fp_xvec2 = gfopen("OUTPUT/wavefield_xvec2_for_dumps.txt","w");
            fp_zvec1 = gfopen("OUTPUT/wavefield_zvec1_for_dumps.txt","w");
            fp_zvec2 = gfopen("OUTPUT/wavefield_zvec2_for_dumps.txt","w");
            for (ix = 0; ix < nx; ix++){
                fprintf(fp_xvec1, "%20.6f\n", xvec1[ix]);
                fprintf(fp_xvec2, "%20.6f\n", xvec2[ix]);
            }
            for (iz = 0; iz < nz; iz++){
                fprintf(fp_zvec1, "%20.6f\n", zvec1[iz]);
                fprintf(fp_zvec2, "%20.6f\n", zvec2[iz]);
            }
        }


        fclose(fp_xvec1); fclose(fp_xvec2);
        fclose(fp_zvec1); fclose(fp_zvec2);

    }

    /* write the velocity components */
    if (imagetype_wavefield_dumps==SNAP_VELOC) {
        if (use_binary_for_wavefield_dumps) {
            /*For Vx*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vx1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Vx[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

            /*For Vz*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vz1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Vz[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

        }
        /* write ASCII version */
        else {
            /*For Vx*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vx1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Vx[iz*nx+ix]);
                fprintf(fp_field,"\n");
            }
            fclose(fp_field);

            /*For Vz*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Vz1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Vz[iz*nx+ix]);
                fprintf(fp_field,"\n");
            }
            fclose(fp_field);

        }
    }

 /* write the velocity components */
    if (imagetype_wavefield_dumps==SNAP_STRESS) {
        if (use_binary_for_wavefield_dumps) {
            /*For Txx1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txx1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Txx[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

            /*For Tzz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Tzz1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Tzz[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

            /*For Txz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txz1.bin",it);
            fp_field = gfopen(fn,"wb");
            for (ix=0; ix<nx; ix++) {
                for (iz=0; iz<nz; iz++)
                    fwrite(&Txz[iz*nx+ix],sizeof(float),1,fp_field);
            }
            fclose(fp_field);

        }
        /* write ASCII version */
        else {
            /*For Txx1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txx1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Txx[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);

            /*For Tzz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Tzz1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Tzz[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);

            /*For Txz1*/
            snprintf(fn, sizeof(fn),"OUTPUT/wavefield%07d_Txz1.txt",it);
            fp_field = gfopen(fn,"w");
            for (iz=0; iz<nz; iz++) {
                for (ix=0; ix<nx; ix++)
                    fprintf(fp_field, "%-20.6e",Txz[iz*nx+ix]);
                fprintf(fp_field, "\n");
            }
            fclose(fp_field);
        }
    }

    return 0;
}

