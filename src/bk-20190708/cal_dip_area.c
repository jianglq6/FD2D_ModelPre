/***************************************************************************
 *
 * This function is used to calculate the dip and area of the layer in grid.
 *
 * Authors: Luqian Jiang
 *          Wei Zhang
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 3/2019: Original version created by Luqian Jiang
 *          2019.03.13: interp1 (the same as matlab ) function
 *
 ***************************************************************************/

#include "global_used.h"
#include "cal_dip_area.h"

/*************************************************************
 * function: getLayer_of_grid
 * Input: Grid coordinates, layer coordinates(clockwise)
 * Output: the coordinate of grid cutting by the layer
 *
 * History: 2019.05.07 Luqian Jiang
 *************************************************************/
void getLayer_of_grid(float *xvec, float *zvec, float *layer_x, float *layer_z,
                      int nx, int nz, int npoints_interfaces,
                      float *layer_of_grid_x, float *layer_of_grid_z, int *npoints_layer)
{
    int idx_g, idx_l, idz_g, idz_l, id_l;
    float *idx = NULL, *idz = NULL, *id_layer_grid_all = NULL;
    int n_idx, n_idz, n_id_all;
    float *id_layer = NULL;
    int i_idx = 0,i_idz = 0;

    idx = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    idz = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    id_layer_grid_all = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    id_layer = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));

    /* get the index of grids in the layer */
    for (idx_g = 0; idx_g < nx; idx_g++) {
        for (idx_l = 0; idx_l < npoints_interfaces-1; idx_l++) {
            if (   (xvec[idx_g]-layer_x[idx_l])*(xvec[idx_g]-layer_x[idx_l+1]) <= 0
                && ( layer_x[idx_l+1]-layer_x[idx_l]) != 0 ) {
                idx[i_idx] = idx_l + (xvec[idx_g]-layer_x[idx_l])/(layer_x[idx_l+1]-layer_x[idx_l]);
                i_idx++;
            }
        }
    }

    for (idz_g = 0; idz_g < nz; idz_g++) {
        for (idz_l = 0; idz_l < npoints_interfaces-1; idz_l++) {
            if (   (zvec[idz_g]-layer_z[idz_l])*(zvec[idz_g]-layer_z[idz_l+1]) <= 0
                && (layer_z[idz_l+1]-layer_z[idz_l]) != 0 ) {
                idz[i_idz] = idz_l + (zvec[idz_g]-layer_z[idz_l])/(layer_z[idz_l+1]-layer_z[idz_l]);
                i_idz++;
            }
        }
    }

    /* sort the index and remove the duplicate items */
    n_idx = i_idx;
    n_idz = i_idz;

    /* Sort and merge */
    shellsort(idx,n_idx);
    shellsort(idz,n_idz);
    merge(idx, idz, n_idx, n_idz, id_layer_grid_all, &n_id_all);


    /* get the coordinate of grids in the layer */
    // TODO !! maybe we can classified these points
    for (id_l = 0; id_l <  npoints_interfaces; id_l++)
        id_layer[id_l] = (float)id_l;

    interp1(id_layer, npoints_interfaces, layer_x, id_layer_grid_all, n_id_all, layer_of_grid_x);
    interp1(id_layer, npoints_interfaces, layer_z, id_layer_grid_all, n_id_all, layer_of_grid_z);

    *npoints_layer = n_id_all;

    free(idx);
    free(idz);
    free(id_layer);
    free(id_layer_grid_all);
}

/******************************************************************
 * calculate the dip and areas of the layer cutting grid
 * !! Can not used in variable mesh
 * !! Must first ensure the layer_of_grid is correct and complete
 * 2019.05.27
 ******************************************************************/
void cal_dip_area(float *xvec, float *zvec, float dx, float dz, int nx, int nz,
                  float *layer_x, float *layer_z, int npoints_interfaces, int *npoints_para,
                  int *ix, int *iz, int *loc_type, float *area1, float *area2, float *theta)
{
    int   ix1, iz1, ix2, iz2, ipoint_layer, npoints_layer,i_point;
    float *layer_of_grid_x = NULL, *layer_of_grid_z = NULL;

    layer_of_grid_x = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    layer_of_grid_z = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));

    /* get the layer of grid */
    getLayer_of_grid(xvec, zvec, layer_x, layer_z, nx, nz,
        npoints_interfaces, layer_of_grid_x, layer_of_grid_z, &npoints_layer);

    *npoints_para = npoints_layer - 1;

    for (i_point = 0; i_point < npoints_layer-1; i_point++) {

        // TODO ?? Precision problem, for example, floor(1.99999...) = 1, may make things wrong!
        ix1 = floor( (layer_of_grid_x[i_point]  -xvec[0]) / dx + SYS_ERR*10 );
        iz1 = floor( (layer_of_grid_z[i_point]  -zvec[0]) / dz + SYS_ERR*10 );

        ix2 = floor( (layer_of_grid_x[i_point+1]-xvec[0]) / dx + SYS_ERR*10 );
        iz2 = floor( (layer_of_grid_z[i_point+1]-zvec[0]) / dz + SYS_ERR*10 );

        printf("i_point: %d, ix1: %d; ix2: %d; iz1: %d; iz2:%d \n",
                i_point, ix1, ix2,iz1, iz2);
        /* Return the index of the upper left point */

        // Check if the interpolation is right!
        if (abs(ix1-ix2) > 1 || abs(iz1-iz2) > 1) {
            printf("layer_of_grid is incorrect! Please check out!\n\n");
            return;
        }

        ix[i_point] = (ix1 <= ix2) ? ix1 : ix2;
        iz[i_point] = (iz1 <= iz2) ? iz1 : iz2;

        printf("xvec[%d]: %f; zvec[%d]: %f; [%f, %f], [%f, %f]\n",
                ix[i_point], xvec[ix[i_point]], iz[i_point], zvec[iz[i_point]],
                layer_of_grid_x[i_point], layer_of_grid_z[i_point],layer_of_grid_x[i_point+1], layer_of_grid_z[i_point+1]);

        theta[i_point] = atan( (layer_of_grid_z[i_point+1] - layer_of_grid_z[i_point]) /
                               (layer_of_grid_x[i_point+1] - layer_of_grid_x[i_point]) );


        if ( abs(layer_of_grid_x[i_point] - xvec[ix[i_point]]) <= SYS_ERR && iz1 == iz[i_point] ) {

            if ( abs(layer_of_grid_z[i_point+1] - zvec[iz[i_point]] <= SYS_ERR && ix2 == ix[i_point]) ) {
                area1[i_point] = ( layer_of_grid_x[i_point+1] - xvec[ix[i_point]] ) *
                                 ( layer_of_grid_z[i_point  ] - zvec[iz[i_point]] )/2;
                area2[i_point] = dx*dz - area1[i_point];
                loc_type[i_point] = 1;
            }
            else if ( abs(layer_of_grid_x[i_point+1] - xvec[ix[i_point]+1]) <= SYS_ERR && iz2 == iz[i_point]) {
                area1[i_point] = ( (layer_of_grid_z[i_point  ] - zvec[iz[i_point]] )
                                  +(layer_of_grid_z[i_point+1] - zvec[iz[i_point]] ) )* dx/2;
                area2[i_point] = dx*dz - area1[i_point];
                loc_type[i_point] = 2;
            }
            else if ( abs(layer_of_grid_z[i_point+1] - zvec[iz[i_point]+1]) <= SYS_ERR && ix2 == ix[i_point]) {
                area2[i_point] = ( layer_of_grid_x[i_point+1] - xvec[ix[i_point]  ] ) *
                                 (-layer_of_grid_z[i_point  ] + zvec[iz[i_point]+1] )/2;
                area1[i_point] = dx*dz - area2[i_point];
                loc_type[i_point] = 3;
            }

        }

        else if ( abs(layer_of_grid_z[i_point] - zvec[iz[i_point]]) <= SYS_ERR && ix1 == ix[i_point] ) {

            if ( abs(layer_of_grid_x[i_point+1] - xvec[ix[i_point]+1]) <= SYS_ERR && iz2 == iz[i_point]) {
                area1[i_point] = (-layer_of_grid_x[i_point  ] + xvec[ix[i_point]+1] ) *
                                 ( layer_of_grid_z[i_point+1] - zvec[iz[i_point]  ] ) /2;
                area2[i_point] = dx*dz - area1[i_point];
                loc_type[i_point] = 4;
            }
            else if ( abs(layer_of_grid_z[i_point+1] - zvec[iz[i_point]+1]) <= SYS_ERR && ix2 == ix[i_point]) {
                area1[i_point] = ( (xvec[ix[i_point]+1] - layer_of_grid_x[i_point  ]  )
                                  +(xvec[ix[i_point]+1] - layer_of_grid_x[i_point+1]  ) ) * dz/2;
                area2[i_point] = dx*dz - area1[i_point];
                loc_type[i_point] = 5;
                printf("%f %f\n",(xvec[ix[i_point]+1] - layer_of_grid_x[i_point  ]  ), (xvec[ix[i_point]+1] - layer_of_grid_x[i_point+1]  ) );
            }
        }

        else if ( abs(layer_of_grid_z[i_point] - zvec[iz[i_point]+1]) <= SYS_ERR && ix1 == ix[i_point] ) {

            if (abs(layer_of_grid_x[i_point+1] - xvec[ix[i_point]+1]) <= SYS_ERR && iz2 == iz[i_point] ) {
                area2[i_point] = ( xvec[ix[i_point]+1] - layer_of_grid_x[i_point  ] ) *
                                 ( zvec[iz[i_point]+1] - layer_of_grid_z[i_point+1] ) / 2;
                area1[i_point] = dx * dz - area2[i_point];
                loc_type[i_point] = 6;
            }

        }

        if ( abs(layer_of_grid_x[i_point+1] - xvec[ix[i_point]]) <= SYS_ERR && iz2 == iz[i_point] ) {

            if ( abs(layer_of_grid_z[i_point] - zvec[iz[i_point]] <= SYS_ERR && ix1 == ix[i_point]) ) {
                area2[i_point] = ( layer_of_grid_x[i_point  ] - xvec[ix[i_point]] ) *
                                 ( layer_of_grid_z[i_point+1] - zvec[iz[i_point]] )/2;
                area1[i_point] = dx*dz - area2[i_point];
                loc_type[i_point] = 7;
            }
            else if ( abs(layer_of_grid_x[i_point] - xvec[ix[i_point]+1]) <= SYS_ERR && iz1 == iz[i_point]) {
                area2[i_point] = ( ( layer_of_grid_z[i_point+1] - zvec[iz[i_point]] )
                                  +( layer_of_grid_z[i_point  ] - zvec[iz[i_point]] ) )* dx/2;
                area1[i_point] = dx*dz - area2[i_point];
                loc_type[i_point] = 8;
            }
            else if ( abs(layer_of_grid_z[i_point] - zvec[iz[i_point]+1]) <= SYS_ERR && ix1 == ix[i_point]) {
                area1[i_point] = ( layer_of_grid_x[i_point  ] - xvec[ix[i_point]  ] ) *
                                 (-layer_of_grid_z[i_point+1] + zvec[iz[i_point]+1] )/2;
                area2[i_point] = dx*dz - area1[i_point];
                loc_type[i_point] = 9;
            }

        }

        else if ( abs(layer_of_grid_z[i_point+1] - zvec[iz[i_point]]) <= SYS_ERR && ix2 == ix[i_point] ) {

            if ( abs(layer_of_grid_x[i_point] - xvec[ix[i_point]+1]) <= SYS_ERR && iz1 == iz[i_point]) {
                area2[i_point] = (-layer_of_grid_x[i_point+1] + xvec[ix[i_point]+1] ) *
                                 ( layer_of_grid_z[i_point  ] - zvec[iz[i_point]  ] ) /2;
                area1[i_point] = dx*dz - area2[i_point];
                loc_type[i_point] = 10;
            }
            else if ( abs(layer_of_grid_z[i_point] - zvec[iz[i_point]+1]) <= SYS_ERR && ix1 == ix[i_point]) {
                area2[i_point] = ( (xvec[ix[i_point]+1] - layer_of_grid_x[i_point  ]  )
                                  +(xvec[ix[i_point]+1] - layer_of_grid_x[i_point+1]  ) ) * dz/2;
                area1[i_point] = dx*dz - area2[i_point];
                loc_type[i_point] = 11;
            }

        }

        else if ( abs(layer_of_grid_z[i_point+1] - zvec[iz[i_point]+1]) <= SYS_ERR && ix2 == ix[i_point] ) {

            if (abs(layer_of_grid_x[i_point] - xvec[ix[i_point]+1]) <= SYS_ERR && iz1 == iz[i_point] ) {
                area1[i_point] = ( xvec[ix[i_point]+1] - layer_of_grid_x[i_point+1] ) *
                                 ( zvec[iz[i_point]+1] - layer_of_grid_z[i_point  ] ) / 2;
                area2[i_point] = dx * dz - area1[i_point];
                loc_type[i_point] = 12;
            }

        }
        printf("area1: %f, area2:%f, loc_type: %d \n",
                area1[i_point],area2[i_point], loc_type[i_point]);
    }

    free(layer_of_grid_x);
    free(layer_of_grid_z);
}


