#include <stdbool.h>
#include "pre_global_used.h"
#include "pre_read_conf.h"
#include "pre_media_parameterization.h"

/* malloc a float array and initialization */
float *creat_para_array(int n)
{
    int i;
    float *array;
    array = (float*) malloc( n * sizeof(float));
    for (i = 0; i < n; i++ ) {
        array[i] = -1.0;
    }
    return array;
}

int main(int agrc, char *argv)
{

    char  config_file[] =  "../DATA/Par_file";

    /* Grid information */
    float xmin, dx, zmin, dz;
    char  xmin_key[] = "xmin", nx_key[] = "nx", dx_key[] = "dx", dz_key[] = "dz", zmin_key[] = "zmin", nz_key[] = "nz";
    int   nx, nz, ix, iz;

    int   nbmodels, effective_para_method;
    char  nbmodels_key[] = "nbmodels", effective_para_method_key[] = "effective_para_method";
    bool  use_existing_model;
    char  use_existing_model_key[] = "use_existing_model";
    char  interfacesfile[MAX_VAL_LEN], interfacesfile_key[] = "interfacesfile";
    char  materialfile[MAX_VAL_LEN], materialfile_key[] = "materialfile";

    float *xvec = NULL, *zvec = NULL;


    /* interfaces */
    int    number_of_interfaces;
    struct interfaces *interface;

    float xmax, zmax;

    /* Read the parameters I need */
    read_conf_value_float(config_file, xmin_key, &xmin);
    read_conf_value_float(config_file, zmin_key, &zmin);
    read_conf_value_float(config_file, dx_key  , &dx  );
    read_conf_value_float(config_file, dz_key  , &dz  );
    read_conf_value_int(config_file, nx_key, &nx);
    read_conf_value_int(config_file, nz_key, &nz);
    read_conf_value_int(config_file, nbmodels_key, &nbmodels);
    read_conf_value_int(config_file, effective_para_method_key, &effective_para_method);
    read_conf_value_bool(config_file, use_existing_model_key, &use_existing_model);

    /******************************************************************
     * Grid Setting
     ******************************************************************/
    xvec = (float*) malloc( nx*sizeof(float) );
    zvec = (float*) malloc( nz*sizeof(float) );
    for (ix = 0; ix<nx; ix++) {
        xvec[ix] = xmin + ix * dx;
    }
    xmax = xvec[nx-1];
    for (iz = 0; iz<nz; iz++) {
        zvec[iz] = zmin + iz * dz;
    }
    zmax = zvec[nz-1];


    /* If not using an existing model, need to be sorted as a function */
    if (!use_existing_model) {
            /* Material */
        int i = 0, j;
        int   *model_number, *material_type, *ierr;

        float *val1 = NULL, *val2 = NULL, *val3 = NULL, *val4 = NULL, *val6 = NULL, *val7 = NULL, *val5 = NULL;
        float *val8 = NULL, *val9 = NULL, *val10 = NULL, *val11 = NULL, *val12 = NULL, *val13 = NULL;
        float *VP = NULL, *VS = NULL, *RHO = NULL, *MU = NULL, *LAM = NULL, *L2M = NULL;      /* model parameter with layer */
        float *lam = NULL, *mu = NULL, *lam2mu = NULL, *muxz = NULL, *Bx = NULL, *Bz = NULL;
        float *rho_x = NULL, *rho_z = NULL;
        float *lam11 = NULL, *mu00 = NULL, *mu11 = NULL, *rho01 = NULL, *rho10 = NULL; /* For LOC, GRI and VOL */
        /* For TTI Approximation */
        float *lam2mu00 = NULL, *lam2mu11 = NULL, *lam00 = NULL;
        float *c11_1 = NULL, *c13_1 = NULL, *c15_1 = NULL, *c33_1 = NULL, *c35_1 = NULL, *c55_1 = NULL;
        float *c11_2 = NULL, *c13_2 = NULL, *c15_2 = NULL, *c33_2 = NULL, *c35_2 = NULL, *c55_2 = NULL;

        FILE *fp1, *fp2, *fp3;

        ierr = &i;

        /* Read material file */
        model_number  = (int*) malloc( nbmodels * sizeof(int) );
        material_type = (int*) malloc( nbmodels * sizeof(int) );
        val1  = (float*) malloc( nbmodels * sizeof(float) );
        val2  = (float*) malloc( nbmodels * sizeof(float) );
        val3  = (float*) malloc( nbmodels * sizeof(float) );
        val4  = (float*) malloc( nbmodels * sizeof(float) );
        val5  = (float*) malloc( nbmodels * sizeof(float) );
        val6  = (float*) malloc( nbmodels * sizeof(float) );
        val7  = (float*) malloc( nbmodels * sizeof(float) );
        val8  = (float*) malloc( nbmodels * sizeof(float) );
        val9  = (float*) malloc( nbmodels * sizeof(float) );
        val10 = (float*) malloc( nbmodels * sizeof(float) );
        val11 = (float*) malloc( nbmodels * sizeof(float) );
        val12 = (float*) malloc( nbmodels * sizeof(float) );
        val13 = (float*) malloc( nbmodels * sizeof(float) );

        /* model parameter (layer) */
        VP  = (float*) malloc( nbmodels * sizeof(float) );
        VS  = (float*) malloc( nbmodels * sizeof(float) );
        RHO = (float*) malloc( nbmodels * sizeof(float) );
        MU  = (float*) malloc( nbmodels * sizeof(float) );
        LAM = (float*) malloc( nbmodels * sizeof(float) );
        L2M = (float*) malloc( nbmodels * sizeof(float) );

        /* For LOC, GRI and VOL */
        lam    = creat_para_array( nx*nz );
        mu     = creat_para_array( nx*nz );
        lam2mu = creat_para_array( nx*nz );
        muxz   = creat_para_array( nx*nz );
        rho_x  = creat_para_array( nx*nz );
        rho_z  = creat_para_array( nx*nz );
        Bx     = creat_para_array( nx*nz );
        Bz     = creat_para_array( nx*nz );

        lam11  = creat_para_array( nx*nz );
        mu00   = creat_para_array( nx*nz );
        mu11   = creat_para_array( nx*nz );
        rho01  = creat_para_array( nx*nz );
        rho10  = creat_para_array( nx*nz );

        c11_1 = creat_para_array( nx*nz );
        c13_1 = creat_para_array( nx*nz );
        c15_1 = creat_para_array( nx*nz );
        c33_1 = creat_para_array( nx*nz );
        c35_1 = creat_para_array( nx*nz );
        c55_1 = creat_para_array( nx*nz );
        c11_2 = creat_para_array( nx*nz );
        c13_2 = creat_para_array( nx*nz );
        c15_2 = creat_para_array( nx*nz );
        c33_2 = creat_para_array( nx*nz );
        c35_2 = creat_para_array( nx*nz );
        c55_2 = creat_para_array( nx*nz );

        lam00    = creat_para_array( nx*nz );
        lam2mu00 = creat_para_array( nx*nz );
        lam2mu11 = creat_para_array( nx*nz );


        /* Initialization */
        read_conf_value_string(config_file, materialfile_key  , materialfile);
        read_conf_value_string(config_file, interfacesfile_key, interfacesfile);

        read_material(nbmodels, materialfile, material_type, model_number,
                val1, val2, val3, val4, val5, val6, val7, val8, val9,
                val10, val11, val12, val13, ierr);

        /* Read interfaces file */
        interface = read_interfaces(interfacesfile, &number_of_interfaces);


        // Need to be modefied!
        for (i = 0; i<nbmodels; i++) {
            switch(material_type[i]) {
                /* elastic media */
                case 1:
                    fprintf(stdout, "material %d is elastic media\n", i+1);
                    RHO[i] = val1[i] ;
                    VP[i]  = val2[i] ;
                    VS[i]  = val3[i] ;
                    MU[i]  = pow(VS[i],2) * RHO[i];
                    L2M[i] = pow(VP[i],2) * RHO[i];
                    LAM[i] = L2M[i] - 2.0 * MU[i];

                    printf("vp: %f, vs: %f, rho: %f, mu: %f, l2mu:%f, lam: %f\n", VP[i], VS[i], RHO[i], MU[i], L2M[i], LAM[i]);
                    break;
            }
        }

        // !! TODO: Maybe we need to do it in the interface file
        /* To avoid that the half point have no value after assignment */

        for (i = 0; i < number_of_interfaces; i++) {
            for (j = 0; j < interface[i].npoints_interfaces; j++) {
                if ( (interface[i].x_loc[j] - xmin) <= 0 && (interface[i].x_loc[j] - xmin) > -dx/2 )
                    interface[i].x_loc[j] = interface[i].x_loc[j] - dx/2;
                if ( (interface[i].x_loc[j] - xmax) >= 0 && (interface[i].x_loc[j] - xmax) < dx/2 )
                    interface[i].x_loc[j] = interface[i].x_loc[j] + dx/2;
                if ( (interface[i].z_loc[j] - zmin) <= 0 && (interface[i].z_loc[j] - zmin) > -dz/2 )
                    interface[i].z_loc[j] = interface[i].z_loc[j] - dz/2;
                if ( (interface[i].z_loc[j] - zmax) >= 0 && (interface[i].z_loc[j] - zmax) < dz/2 )
                    interface[i].z_loc[j] = interface[i].z_loc[j] + dz/2;

                printf("%f %f\n", interface[i].x_loc[j], interface[i].z_loc[j]);
            }
            printf("\n");
        }

        /* Elastic media: Assignment && Parameterization && Forward */
        switch(effective_para_method) {

            case 0:
                fprintf(stdout, "LOC: DO NOT use effective media method.\n");
                media_parameterization_loc(interface, number_of_interfaces,
                    xmin, zmin, dx, dz, nx, nz, xvec, zvec,
                    LAM, MU, RHO, L2M,
                    lam, mu, muxz, rho_x, rho_z);

                fp1 = fopen("lam_loc.dat","w");
                for (i = 0; i < nx; i++)
                    for (j = 0; j< nz; j++) {
                        fwrite(&lam[j*nx+i],sizeof(float),1,fp1);
                    }
                fclose(fp1);

                fp2 = fopen("mu_loc.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j< nz; j++) {
                        fwrite(&mu[j*nx+i],sizeof(float),1,fp2);
                        //printf("%f ",mu[j*nx+i]);
                    }
                }
                fclose(fp2);

                fp3 = fopen("muxz_loc.dat","w");
                for (i = 0; i < nx; i++)
                    for (j = 0; j< nz; j++) {
                        fwrite(&muxz[j*nx+i],sizeof(float),1,fp3);
                    }
                fclose(fp3);

                break;

            case 1:
                fprintf(stdout, "GRI: Grid average method.\n");
                media_parameterization_gri(interface, number_of_interfaces,
                    xmin, zmin, dx, dz, nx, nz, xvec, zvec,
                    LAM, MU, RHO, L2M,
                    lam, lam2mu, muxz, Bx, Bz, rho_x, mu);
                break;

            case 2:
                fprintf(stdout, "VOL: volume arithmetic and harmonic averaging method.\n");
                media_parameterization_vol(interface, number_of_interfaces,
                    xmin, zmin, dx, dz, nx, nz, xvec, zvec,
                    LAM, MU, RHO, lam11, mu00, mu11, rho01, rho10,
                    lam, mu, muxz, rho_x, rho_z);

                fp1 = fopen("lam_vol.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&lam[j*nx+i],sizeof(float),1,fp1);
                    }
                }
                fclose(fp1);


                fp2 = fopen("mu_vol.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&mu[j*nx+i],sizeof(float),1,fp2);
                        //printf("%f  ",mu[j*nx+i]);
                    }
                    //printf("\n");
                }
                fclose(fp2);

                fp3 = fopen("muxz_vol.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&muxz[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);


                fp3 = fopen("rhox_vol.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&rho_x[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("rhoz_vol.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&rho_z[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                break;
            case 4:
                fprintf(stdout, "TTI: TTI approximation method.\n");
                media_parameterization_tti(interface, number_of_interfaces,
                     xmin, zmin, dx, dz, nx, nz, xvec, zvec,
                     LAM, MU, RHO, L2M,
                     lam2mu00, lam2mu11, lam00, lam11, mu00, mu11, rho01, rho10,
                     c11_1, c13_1, c15_1, c33_1, c35_1, c55_1,
                     c11_2, c13_2, c15_2, c33_2, c35_2, c55_2,
                     rho_x, rho_z);

                fp1 = fopen("c11_1.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c11_1[j*nx+i],sizeof(float),1,fp1);
                    }
                }
                fclose(fp1);


                fp2 = fopen("c55_1.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c55_1[j*nx+i],sizeof(float),1,fp2);
                        //printf("%f  ",mu[j*nx+i]);
                    }
                    //printf("\n");
                }
                fclose(fp2);

                fp3 = fopen("c13_1.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c13_1[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("c35_1.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c35_1[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("c33_1.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c33_1[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("c15_1.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c15_1[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp1 = fopen("c11_2.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c11_2[j*nx+i],sizeof(float),1,fp1);
                    }
                }
                fclose(fp1);


                fp2 = fopen("c55_2.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c55_2[j*nx+i],sizeof(float),1,fp2);
                        //printf("%f  ",mu[j*nx+i]);
                    }
                    //printf("\n");
                }
                fclose(fp2);

                fp3 = fopen("c13_2.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c13_2[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("c35_2.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c35_2[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("c33_2.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c33_2[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                fp3 = fopen("c15_2.dat","w");
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nz; j++) {
                        fwrite(&c15_2[j*nx+i],sizeof(float),1,fp3);
                    }
                }
                fclose(fp3);

                break;

            default:
                fprintf(stdout,"You should fill 0-4");
                break;

            free(model_number);
            free(material_type);
            free(val1);
            free(val2);
            free(val3);
            free(val4);
            free(val5);
            free(val6);
            free(val7);
            free(val8);
            free(val9);
            free(val10);
            free(val11);
            free(val12);
            free(val13);

            free(VP);
            free(VS);
            free(RHO);
            free(MU);
            free(LAM);
            free(L2M);

            free(lam);
            free(mu);
            free(lam2mu);
            free(muxz);
            free(rho_x);
            free(rho_z);
            free(Bx);
            free(Bz);

            free(lam11);
            free(mu00);
            free(mu11);
            free(rho01);
            free(rho10);

            free(c11_1);
            free(c13_1);
            free(c15_1);
            free(c35_1);
            free(c33_1);
            free(c55_1);
            free(c11_2);
            free(c13_2);
            free(c15_2);
            free(c35_2);
            free(c33_2);
            free(c55_2);
            free(lam2mu00);
            free(lam2mu11);
            free(lam00);

        }

    }


    free(xvec);
    free(zvec);

    return 0;
}

