/***************************************************************************
 *
 * This function is used to calculate the dip and area of the layer in grid.
 *
 * Authors: Luqian Jiang
 *          Wei Zhang
 *
 * Copyright (c) 2018 - 2019 zwlab
 *
 * History: 03/2019: Original version created by Luqian Jiang
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
    int   ix1, iz1, ix2, iz2, ipoint_layer, npoints_layer, i_point, i_point_para = 0;
    float *layer_of_grid_x = NULL, *layer_of_grid_z = NULL;
    int ix_para, iz_para, loc_type_para;
    float area1_para, area2_para, theta_para;
    FILE *fp;
    fp = fopen("layer_grid.dat","w");

    layer_of_grid_x = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));
    layer_of_grid_z = (float*)malloc(MAX_LAYER_POINTS*sizeof(float));

    /* get the layer of grid */
    getLayer_of_grid(xvec, zvec, layer_x, layer_z, nx, nz,
        npoints_interfaces, layer_of_grid_x, layer_of_grid_z, &npoints_layer);

    //*npoints_para = npoints_layer - 1;

    for (i_point = 0; i_point < npoints_layer-1; i_point++) {

        // TODO ?? Precision problem, for example, floor(1.99999...) = 1, may make things wrong!
        ix1 = floor( (layer_of_grid_x[i_point]  -xvec[0]) / dx + 1e-3 );
        iz1 = floor( (layer_of_grid_z[i_point]  -zvec[0]) / dz + 1e-3 );

        ix2 = floor( (layer_of_grid_x[i_point+1]-xvec[0]) / dx + 1e-3 );
        iz2 = floor( (layer_of_grid_z[i_point+1]-zvec[0]) / dz + 1e-3 );

        /* Return the index of the upper left point */

        // Check if the interpolation is right!
        if (abs(ix1-ix2) > 1 || abs(iz1-iz2) > 1 && iz1 ) {
            printf("layer_of_grid is incorrect! Please check out!\n");
            printf("ix1: %d, ix2: %d, iz1: %d, iz2: %d, layer_of_grid_x[%d]: (%f %f), layer_of_grid_z[%d]: (%f %f)\n",
                    ix1,ix2,iz1,iz2,i_point,layer_of_grid_x[i_point],layer_of_grid_z[i_point], i_point+1, layer_of_grid_x[i_point+1], layer_of_grid_z[i_point+1]);
            //return;
        }

        ix_para = (ix1 <= ix2) ? ix1 : ix2;
        iz_para = (iz1 <= iz2) ? iz1 : iz2;


        theta_para = atan2( (layer_of_grid_z[i_point+1] - layer_of_grid_z[i_point]),
                            (layer_of_grid_x[i_point+1] - layer_of_grid_x[i_point]) );


        /* 12 types of the layer pass through the grid */
        if ( fabs(layer_of_grid_x[i_point  ]-xvec[ix_para]) <= 8e-4*dx && ( iz1 == iz_para || iz1 == iz_para+1 ) &&
             fabs(layer_of_grid_z[i_point+1]-zvec[iz_para]) <= 8e-4*dz && ( ix2 == ix_para || ix2 == ix_para+1 ) ){
            area1_para = ( layer_of_grid_x[i_point+1] - xvec[ix_para] ) *
                         ( layer_of_grid_z[i_point  ] - zvec[iz_para] )/2;
            area2_para = dx*dz - area1_para;
            loc_type_para = 1;
        }
        else if ( fabs(layer_of_grid_z[i_point  ]-zvec[iz_para  ]) <= 8e-4*dz && (ix1 == ix_para || ix1 == ix_para+1) &&
                  fabs(layer_of_grid_x[i_point+1]-xvec[ix_para+1]) <= 8e-4*dx && (iz2 == iz_para || iz2 == iz_para+1) ) {
            area1_para = (-layer_of_grid_x[i_point  ] + xvec[ix_para+1] ) *
                         ( layer_of_grid_z[i_point+1] - zvec[iz_para  ] ) /2;
            area2_para = dx*dz - area1_para;
            loc_type_para = 4;
        }
        else if ( fabs(layer_of_grid_x[i_point  ]-xvec[ix_para  ]) <= 8e-4*dx && ( iz1 == iz_para || iz1 == iz_para+1 ) &&
                  fabs(layer_of_grid_x[i_point+1]-xvec[ix_para+1]) <= 8e-4*dx && ( iz2 == iz_para || iz2 == iz_para+1 ) ) {
            area1_para = ( (layer_of_grid_z[i_point  ] - zvec[iz_para] )
                          +(layer_of_grid_z[i_point+1] - zvec[iz_para] ) )* dx/2;
            area2_para = dx*dz - area1_para;
            loc_type_para = 2;
        }
        else if ( fabs(layer_of_grid_x[i_point  ]-xvec[ix_para  ]) <= 8e-4*dx && (iz1 == iz_para || iz1 == iz_para+1) &&
                  fabs(layer_of_grid_z[i_point+1]-zvec[iz_para+1]) <= 8e-4*dz && (ix2 == ix_para ) ) {
            area2_para = ( layer_of_grid_x[i_point+1] - xvec[ix_para  ] ) *
                         (-layer_of_grid_z[i_point  ] + zvec[iz_para+1] )/2;
            area1_para = dx*dz - area2_para;
            loc_type_para = 3;
        }
        else if ( fabs(layer_of_grid_x[i_point  ]-xvec[ix_para+1]) <= 8e-4*dx && (iz1 == iz_para || iz1 == iz_para+1) &&
                  fabs(layer_of_grid_z[i_point+1]-zvec[iz_para+1]) <= 8e-4*dz && (ix2 == ix_para || ix2 == ix_para+1) ) {
            area1_para = ( xvec[ix_para+1] - layer_of_grid_x[i_point+1] ) *
                         ( zvec[iz_para+1] - layer_of_grid_z[i_point  ] ) / 2;
            area2_para = dx * dz - area1_para;
            loc_type_para = 6;
        }
        else if ( fabs(layer_of_grid_z[i_point  ]-zvec[iz_para  ]) <= 8e-4*dz && (ix1 == ix_para )  &&
                  fabs(layer_of_grid_z[i_point+1]-zvec[iz_para+1]) <= 8e-4*dz && (ix2 == ix_para ) ) {
             area1_para = ( (xvec[ix_para+1] - layer_of_grid_x[i_point  ]  )
                           +(xvec[ix_para+1] - layer_of_grid_x[i_point+1]  ) ) * dz/2;
             area2_para = dx*dz - area1_para;
             loc_type_para = 5;
        }
        else if ( fabs(layer_of_grid_x[i_point+1]-xvec[ix_para]) <= 8e-4*dx && (iz2 == iz_para ) &&
                  fabs(layer_of_grid_z[i_point  ]-zvec[iz_para]) <= 8e-4*dz && (ix1 == ix_para ) ) {
            area2_para = ( layer_of_grid_x[i_point  ] - xvec[ix_para] ) *
                         ( layer_of_grid_z[i_point+1] - zvec[iz_para] )/2;
            area1_para = dx*dz - area2_para;
            loc_type_para = 7;
        }
        else if ( fabs(layer_of_grid_x[i_point+1]-xvec[ix_para  ]) <= 8e-4*dx && (iz2 == iz_para ) &&
                  fabs(layer_of_grid_x[i_point  ]-xvec[ix_para+1]) <= 8e-4*dx && (iz1 == iz_para ) ) {
            area2_para = ( ( layer_of_grid_z[i_point+1] - zvec[iz_para] )
                          +( layer_of_grid_z[i_point  ] - zvec[iz_para] ) )* dx/2;
            area1_para = dx*dz - area2_para;
            loc_type_para = 8;
        }
        else if ( fabs(layer_of_grid_x[i_point+1]-xvec[ix_para  ]) <= 8e-4*dx && (iz2 == iz_para || iz2 == iz_para+1) &&
                  fabs(layer_of_grid_z[i_point  ]-zvec[iz_para+1]) <= 8e-4*dz && (ix1 == ix_para || ix1 == ix_para+1) ) {
            area1_para = ( layer_of_grid_x[i_point  ] - xvec[ix_para  ] ) *
                         (-layer_of_grid_z[i_point+1] + zvec[iz_para+1] )/2;
            area2_para = dx*dz - area1_para;
            loc_type_para = 9;
        }
        else if ( fabs(layer_of_grid_z[i_point+1]-zvec[iz_para  ]) <= 8e-4*dz && (ix2 == ix_para || ix2 == ix_para+1) &&
                  fabs(layer_of_grid_x[i_point  ]-xvec[ix_para+1]) <= 8e-4*dx && (iz1 == iz_para                        ) ) {
            area2_para = (-layer_of_grid_x[i_point+1] + xvec[ix_para+1] ) *
                         ( layer_of_grid_z[i_point  ] - zvec[iz_para  ] ) /2;
            area1_para = dx*dz - area2_para;
            loc_type_para = 10;
        }
        else if ( fabs(layer_of_grid_z[i_point+1]-zvec[iz_para  ]) <= 8e-4*dz && (ix2 == ix_para || ix2 == ix_para+1)  &&
                  fabs(layer_of_grid_z[i_point  ]-zvec[iz_para+1]) <= 8e-4*dz && (ix1 == ix_para || ix1 == ix_para+1) ) {
            area2_para = ( (xvec[ix_para+1] - layer_of_grid_x[i_point  ]  )
                          +(xvec[ix_para+1] - layer_of_grid_x[i_point+1]  ) ) * dz/2;
            area1_para = dx*dz - area2_para;
            loc_type_para = 11;
        }
        else if ( fabs(layer_of_grid_x[i_point+1]-xvec[ix_para+1]) <= 8e-4*dx && (iz2 == iz_para || iz2 == iz_para+1) &&
                  fabs(layer_of_grid_z[i_point  ]-zvec[iz_para+1]) <= 8e-4*dz && (ix1 == ix_para || ix1 == ix_para+1) ) {
            area2_para = ( xvec[ix_para+1] - layer_of_grid_x[i_point  ] ) *
                         ( zvec[iz_para+1] - layer_of_grid_z[i_point+1] ) / 2;
            area1_para = dx * dz - area2_para;
            loc_type_para = 12;
        }
        else {
            loc_type_para = 0;
            if (iz_para > 0 && iz_para < nz-1 && ix_para > 0 && ix_para < nx-1 ) {
                if ( fabs(layer_of_grid_x[i_point]-layer_of_grid_x[i_point+1])>8e-4 &&
                     fabs(layer_of_grid_z[i_point]-layer_of_grid_z[i_point+1])>8e-4)
                    printf("xvec[%d](%f),zvec[%d](%f): The layer_of_grid(%d) has no type, please check!: (%f %f), (%f %f) ix1:%d, ix2: %d, iz1: %d, iz2: %d.   %d\n", ix_para, xvec[ix_para],iz_para,zvec[iz_para],i_point, layer_of_grid_x[i_point],layer_of_grid_z[i_point],layer_of_grid_x[i_point+1], layer_of_grid_z[i_point+1],ix1,ix2,iz1,iz2,loc_type_para);
            }
        }

        /* remove the duplicate point*/
        if (loc_type_para != 0){
            ix[i_point_para] = ix_para;
            iz[i_point_para] = iz_para;
            loc_type[i_point_para] = loc_type_para;
            area1[i_point_para]    = area1_para;
            area2[i_point_para]    = area2_para;
            theta[i_point_para]    = theta_para;
            layer_of_grid_x[i_point_para] = layer_of_grid_x[i_point];
            layer_of_grid_z[i_point_para] = layer_of_grid_z[i_point];
            layer_of_grid_x[i_point_para+1] = layer_of_grid_x[i_point+1];
            layer_of_grid_z[i_point_para+1] = layer_of_grid_z[i_point+1];
            i_point_para ++;
        }

        //if (loc_type_para != 0)
            //fprintf(fp, "%f     %f     %d     %f   %f     %d  %d \n  ",layer_of_grid_x[i_point],layer_of_grid_z[i_point], loc_type_para, area1_para, area2_para, ix_para, iz_para );


    }

    *npoints_para = i_point_para;

    /* Some special circumstances */
    for (i_point = 0; i_point < *npoints_para-1; i_point++) {
        if ( ix[i_point] == ix[i_point+1] && iz[i_point] == iz[i_point+1] ) {
            /* Tangent to the upper side of the grid */
            if ( loc_type[i_point] == 1 && loc_type[i_point+1] == 4 &&
                 layer_of_grid_x[i_point+1] < xvec[ix[i_point]+1] && layer_of_grid_x[i_point+1] > xvec[ix[i_point]] ) {
                loc_type[i_point+1] = 2;
                area1[i_point+1] = area1[i_point] + area1[i_point+1];
                area2[i_point+1] = dx*dz - area1[i_point+1];
                theta[i_point+1] = atan2( (layer_of_grid_z[i_point+2] - layer_of_grid_z[i_point]),
                                          (layer_of_grid_x[i_point+2] - layer_of_grid_x[i_point]) );
            }
            /* Tangent to the right side of the grid */
            if ( loc_type[i_point] == 4 && loc_type[i_point+1] == 6 &&
                 layer_of_grid_z[i_point+1] < zvec[iz[i_point]+1] && layer_of_grid_z[i_point+1] > zvec[iz[i_point]] ) {
                loc_type[i_point+1] = 5;
                area1[i_point+1] = area1[i_point] + area1[i_point+1];
                area2[i_point+1] = dx * dz - area1[i_point+1];
                theta[i_point+1] = atan2( (layer_of_grid_z[i_point+2] - layer_of_grid_z[i_point]),
                                          (layer_of_grid_x[i_point+2] - layer_of_grid_x[i_point]) );
            }
            /* Tangent to the bottom side of the grid */
            if ( loc_type[i_point] == 6 && loc_type[i_point+1] == 9 &&
                 layer_of_grid_x[i_point+1] < xvec[ix[i_point]+1] && layer_of_grid_x[i_point+1] > xvec[ix[i_point]] ) {
                loc_type[i_point+1] = 8;
                area1[i_point+1] = area1[i_point] + area1[i_point+1];
                area2[i_point+1] = dx * dz - area1[i_point+1];
                theta[i_point+1] = atan2( (layer_of_grid_z[i_point+2] - layer_of_grid_z[i_point]),
                                          (layer_of_grid_x[i_point+2] - layer_of_grid_x[i_point]) );
            }
            /* Tangent to the left side of the grid */
            if ( loc_type[i_point] == 9 && loc_type[i_point+1] == 1 &&
                 layer_of_grid_z[i_point+1] < zvec[iz[i_point]+1] && layer_of_grid_z[i_point+1] > zvec[iz[i_point]] ) {
                loc_type[i_point+1] = 11;
                area1[i_point+1] = area1[i_point] + area1[i_point+1];
                area2[i_point+1] = dx * dz - area1[i_point+1];
                theta[i_point+1] = atan( (layer_of_grid_z[i_point+2] - layer_of_grid_z[i_point])/
                                          (layer_of_grid_x[i_point+2] - layer_of_grid_x[i_point]) );
            }
        }



        fprintf(fp, " %f    %f    %d    %f   %f   %d     %d\n  ",layer_of_grid_x[i_point], layer_of_grid_z[i_point],  loc_type[i_point], area1[i_point], area2[i_point], ix[i_point], iz[i_point]);
//
    }



    fclose(fp);

    free(layer_of_grid_x);
    free(layer_of_grid_z);
}



/* Shell sort */
void shellsort(float *arr, int N)
{
    int i,j,Increment;
    float tmp;

    for (Increment = N/2; Increment > 0; Increment /= 2) {
        for (i = Increment; i < N; i++) {
            tmp = arr[i];
            for (j = i; j >= Increment; j -= Increment) {
                if(tmp < arr[j - Increment])
                    arr[j] = arr[j - Increment];
                else
                    break;
            }
            arr[j] = tmp;
        }
    }
}

/* Merge two sorted array and Deduplication (no deduplication -.-) */
void merge(float *arr1, float *arr2, int n1, int n2, float *arr3, int *n3)
{
    int i = 0, j = 0, k = 0;

    // Traverse both array
    while (i < n1 && j < n2) {

        if ( fabs(arr1[i] - arr2[j]) <= 1e-4) {
            arr3[k] = arr2[j];
            k++;
            j++;
            i++;
        }
        else if (arr1[i] < arr2[j]) {
            arr3[k] = arr1[i];
            k++;
            i++;
        }
        else {    // arr1[i] > arr2[j]
            arr3[k] = arr2[j];
            k++;
            j++;
        }

    }

    while (i < n1) {
        arr3[k] = arr1[i];
        k++;
        i++;
    }

    while (j < n2) {
        arr3[k] = arr2[j];
        k++;
        j++;
    }

    *n3 = k;
}

int findNearestNeighborIndex( float value, float *x, int len)
{
    float dist, newDist;
    int idx, i;

    idx = -1;    // expolation: INF(idx = -1) or
    dist = DBL_MAX;
    for(i = 0; i < len; i++) {
        newDist = value - x[i];
        if (newDist >= 0 && newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

// Lagrange interpolation
void interp1(float *x, int x_len, float *v, float *xq, int xq_len, float *vq)
{
    float dx, dy, *slope, *intercept;
    int i, indiceEnVector;

    slope     = (float*) malloc( x_len*sizeof(float) );
    intercept = (float*) malloc( x_len*sizeof(float) );

    for (i = 0; i < x_len; i++) {
        if (i < x_len-1) {
            dx = x[i+1] - x[i];
            dy = v[i+1] - v[i];
            slope[i] = dy / dx;
            intercept[i] = v[i] - x[i] * slope[i];
        } else {
            slope[i] = slope[i-1];
            intercept[i] = intercept[i-1];
        }
    }

    for (i = 0; i < xq_len; i++) {
        indiceEnVector = findNearestNeighborIndex(xq[i], x, x_len);
        //indiceEnVector = (int) xq[i] ;  // For our problem
        if (indiceEnVector != -1)
            vq[i] = slope[indiceEnVector] * xq[i] + intercept[indiceEnVector];
        else
            vq[i] = DBL_MAX;
    }

    free(slope);
    free(intercept);

}

