/*************************************************************************
 *
 * This function is used for model prepare when no use existing model.
 *
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018-2019 zwlab
 *
 * Version: 1.0
 *
 * Date: 11/2019
 *
 * History:
 *     11/2019: Original version created by Luqian Jiang
 *
 ************************************************************************/

#include "pre_interface_struct.h"
#include "read_interface_material.h"
#include "pre_model_prepare.h"
#include "pre_media_parameterization.h"

int model_prepare(int nbmodels, char *materialfile, char *interfacesfile,
        int effective_para_method, int save_model, int *prepared_media,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2,
        float *xvec1, float *zvec1, float *xvec2, float *zvec2, int nx, int nz, float dx, float dz,
        float *c11_1, float *c13_1, float *c15_1, float *c33_1, float *c35_1, float *c55_1,
        float *c11_2, float *c13_2, float *c15_2, float *c33_2, float *c35_2, float *c55_2,
        float *B01, float *B10,
        float *lam2mu00, float *lam2mu11, float *lam00, float *lam11, float *mu00, float *mu11,
        float *rho01, float *rho10, float *rho_x, float *rho_z)
{
    int i, j, ix, iz, ierr = 0;
    int *model_number, *material_type;
    float *val1 = NULL, *val2 = NULL, *val3 = NULL, *val4 = NULL, *val6 = NULL, *val7 = NULL, *val5 = NULL;
    float *val8 = NULL, *val9 = NULL, *val10 = NULL, *val11 = NULL, *val12 = NULL, *val13 = NULL;
    float *VP = NULL, *VS = NULL, *RHO = NULL, *MU = NULL, *LAM = NULL, *L2M = NULL;  /* model parameter with layer */
    float xmin, xmax, zmin, zmax;
    int number_of_interfaces;
    struct interfaces *interface;


    /* For reading material file */
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

    /* read material file */
    read_material(nbmodels, materialfile, material_type, model_number,
            val1, val2, val3, val4, val5, val6, val7, val8, val9,
            val10, val11, val12, val13, &ierr);

    /* Read interfaces file */
    interface = read_interfaces(interfacesfile, &number_of_interfaces);

    // TODO: Not include every situation, just for isotropic media!
    for (i = 0; i<nbmodels; i++) {
        switch(material_type[i]) {
            /* isotropic media */
            case 1:
                fprintf(stdout, "material %d is elastic media\n", i+1);
                RHO[i] = val1[i] ;
                VP[i]  = val2[i] ;
                VS[i]  = val3[i] ;
                MU[i]  = pow(VS[i],2) * RHO[i];
                L2M[i] = pow(VP[i],2) * RHO[i];
                LAM[i] = L2M[i] - 2.0 * MU[i];
                break;

        }
    }

    /* To avoid that the half point have no value after assignment */
    xmin = xvec1[nghost_x1]; xmax = xvec1[nx-nghost_x2-1];
    zmin = zvec1[nghost_z1]; zmax = zvec1[nz-nghost_z2-1];
   for (i = 0; i < number_of_interfaces; i++) {
       for (j = 0; j < interface[i].npoints_interfaces; j++) {
           if ( (interface[i].x_loc[j] - xmin) <= 0 && (interface[i].x_loc[j] - xmin) > -dx/2 )
               interface[i].x_loc[j] = interface[i].x_loc[j] - dx;
           if ( (interface[i].x_loc[j] - xmax) >= 0 && (interface[i].x_loc[j] - xmax) < dx/2 )
               interface[i].x_loc[j] = interface[i].x_loc[j] + dx;
           if ( (interface[i].z_loc[j] - zmin) <= 0 && (interface[i].z_loc[j] - zmin) > -dz/2 )
               interface[i].z_loc[j] = interface[i].z_loc[j] - dz;
           if ( (interface[i].z_loc[j] - zmax) >= 0 && (interface[i].z_loc[j] - zmax) < dz/2 )
               interface[i].z_loc[j] = interface[i].z_loc[j] + dz;
//           printf("%f %f\n", interface[i].x_loc[j], interface[i].z_loc[j]);
       }
//       printf("\n");
   }

    /* effective media parameterization */
    switch(effective_para_method) {

        case LOC:
            fprintf(stdout, "LOC: DO NOT use effective media method.\n");
            ierr = media_parameterization_loc(interface, number_of_interfaces,
                    dx, dz, nx, nz, xvec1, zvec1, xvec2, zvec2,
                    LAM, MU, RHO, L2M,
                    c11_1, c13_1, c55_1, rho_x, rho_z,
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
            /* get B01, B10 */
            for (ix = nghost_x1; ix < nx-nghost_x2; ix ++) {
                for (iz = nghost_z1; iz < nz-nghost_z2; iz++) {
                    B01[iz*nx+ix] = 1/rho_x[iz*nx+ix];
                    B10[iz*nx+ix] = 1/rho_z[iz*nx+ix];
                }
            }
            *prepared_media = PREPARE_ISO;
            break;

        case GRI:
            fprintf(stdout, "GRI: Grid average method.\n");
            ierr = media_parameterization_gri(interface, number_of_interfaces,
                    dx, dz, nx, nz, xvec1, zvec1, xvec2, zvec2,
                    LAM, MU, RHO, L2M,
                    c11_1, c13_1, c55_1, B01, B10, rho_x, mu11,
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
            *prepared_media = PREPARE_ISO;

            break;

        case VOL:
            fprintf(stdout, "VOL: volume arithmetic and harmonic averaging method.\n");
            ierr = media_parameterization_vol(interface, number_of_interfaces,
                    dx, dz, nx, nz, xvec1, zvec1, xvec2, zvec2,
                    LAM, MU, RHO, L2M,
                    lam2mu11, lam11, mu00, rho01, rho10,
                    c11_1, c13_1, c55_1, rho_x, rho_z,
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
            /* get B01, B10 */
            for (ix = nghost_x1; ix < nx-nghost_x2; ix ++) {
                for (iz = nghost_z1; iz < nz-nghost_z2; iz++) {
                    B01[iz*nx+ix] = 1/rho_x[iz*nx+ix];
                    B10[iz*nx+ix] = 1/rho_z[iz*nx+ix];
                }
            }
            *prepared_media = PREPARE_ISO;

            break;

        case TTI:
            fprintf(stdout, "TTI: TTI approximation method.\n");
            ierr = media_parameterization_tti(interface, number_of_interfaces,
                    dx, dz, nx, nz, xvec1, zvec1, xvec2, zvec2,
                    LAM, MU, RHO, L2M,
                    lam2mu00, lam2mu11, lam00, lam11, mu00, mu11, rho01, rho10,
                    c11_1, c13_1, c15_1, c33_1, c35_1, c55_1,
                    c11_2, c13_2, c15_2, c33_2, c35_2, c55_2,
                    rho_x, rho_z,
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
            /* get B01, B10 */
            for (ix = nghost_x1; ix < nx-nghost_x2; ix ++) {
                for (iz = nghost_z1; iz < nz-nghost_z2; iz++) {
                    B01[iz*nx+ix] = 1/rho_x[iz*nx+ix];
                    B10[iz*nx+ix] = 1/rho_z[iz*nx+ix];
                }
            }
            *prepared_media = PREPARE_TTI;

            break;

    }


    /* boundary expansion */
    if (*prepared_media == PREPARE_ISO) {
        ierr = boundary_expansion(c11_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c13_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c55_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(B01, nx, nz,   \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(B10, nx, nz,   \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    }
    else if (*prepared_media == PREPARE_TTI) {
        ierr = boundary_expansion(c11_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c13_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c15_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c33_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c35_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c55_1, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c11_2, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c13_2, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c15_2, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c33_2, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c35_2, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(c55_2, nx, nz, \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(B01, nx, nz,   \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        ierr = boundary_expansion(B10, nx, nz,   \
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
    }


    /* save model */
    if (save_model) {
        if (*prepared_media == PREPARE_ISO) {
            ierr = write_model_iso(save_model, c11_1, c13_1, c55_1, rho_x, rho_z,
                    nx, nz, nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        }
        else if (*prepared_media == PREPARE_TTI) {
            ierr = write_model_TTI(save_model,
                    c11_1, c13_1, c33_1,
                    c15_1, c35_1, c55_1,
                    c11_2, c13_2, c33_2,
                    c15_2, c35_2, c55_2,
                    B01, B10, nx, nz,
                    nghost_x1, nghost_x2, nghost_z1, nghost_z2);
        }
    }



    free(model_number); free(material_type);
    free(val1); free(val2); free(val3); free(val4);
    free(val5); free(val6); free(val7); free(val8);
    free(val9); free(val10); free(val11);
    free(val12); free(val13);

    free(VP); free(VS); free(RHO);
    free(MU); free(LAM); free(L2M);

    return ierr;
}

int boundary_expansion(float *u, int nx, int nz,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    int ix, iz;

    /* top */
    for (ix = 0; ix < nx; ix ++) {
        for (iz = 0; iz < nghost_z1; iz++) {
            u[iz*nx+ix] = u[nghost_z1*nx+ix];
        }
    }


    /* bottom */
    for (ix = 0; ix < nx; ix ++) {
        for (iz = nz-nghost_z2; iz < nz; iz++) {
            u[iz*nx+ix] = u[(nz-nghost_z2-1)*nx+ix];
        }
    }

    /* left */
    for (iz = 0; iz < nz; iz++) {
        for (ix = 0; ix < nghost_x1+1; ix ++) {
            u[iz*nx+ix] = u[iz*nx+nghost_x1+1];
        }
    }

    /* right */
    for (ix = nx-nghost_x2-1; ix < nx; ix ++) {
        for (iz = 0; iz < nz; iz++) {
            u[iz*nx+ix] = u[iz*nx+(nx-nghost_x2-8)];
        }
    }

    return 0;
}


int write_model_iso(int save_model, float *lam2mu, float *lam, float *muxz,
        float *rho_x, float *rho_z, int nx_all, int nz_all,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    int i, j, nx, nz, indx;
    FILE *fp;

    /* save model as binary format */
    if (save_model == 1){
        fp = fopen("OUTPUT/lam.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&lam[j*nx_all+i],sizeof(float),1,fp);
            }
        fclose(fp);
        fp = fopen("OUTPUT/lam2mu.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&lam2mu[j*nx_all+i],sizeof(float),1,fp);
            }
        fclose(fp);
        fp = fopen("OUTPUT/muxz.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&muxz[j*nx_all+i],sizeof(float),1,fp);
            }
        fclose(fp);
        fp = fopen("OUTPUT/rhox.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&rho_x[j*nx_all+i],sizeof(float),1,fp);
            }
        fclose(fp);
        fp = fopen("OUTPUT/rhoz.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&rho_z[j*nx_all+i],sizeof(float),1,fp);
            }
        fclose(fp);
    }

    return 0;
}

int write_model_TTI(int save_model,
        float *c11_1, float *c13_1, float *c33_1,
        float *c15_1, float *c35_1, float *c55_1,
        float *c11_2, float *c13_2, float *c33_2,
        float *c15_2, float *c35_2, float *c55_2,
        float *B01, float *B10, int nx_all, int nz_all,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    int i, j, nx;
    FILE *fp;

    nx = nx_all;// - nghost_x1 - nghost_x2;
    printf("save TTI: %d %d\n", nx_all, nz_all);
    /* save model as binary format */
    if (save_model == 1){

        fp = fopen("OUTPUT/c11_1.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&c11_1[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c13_1.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c13_1[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c55_1.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&c55_1[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c15_1.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c15_1[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c33_1.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c33_1[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c35_1.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c35_1[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c11_2.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&c11_2[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c13_2.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c13_2[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c55_2.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j< nz_all; j++) {
                fwrite(&c55_2[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c15_2.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c15_2[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c33_2.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c33_2[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/c35_2.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&c35_2[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);

        fp = fopen("OUTPUT/B01.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&B01[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);
        fp = fopen("OUTPUT/B10.dat","w");
        for (i = 0; i < nx_all; i++)
            for (j = 0; j < nz_all; j++) {
                fwrite(&B10[j*nx+i],sizeof(float),1,fp);
            }
        fclose(fp);
    }

    return 0;
}
