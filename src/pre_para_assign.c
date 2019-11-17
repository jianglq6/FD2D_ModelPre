/*************************************************************************
 *
 * This function is used to fill the layer closuer line with the media
 * parameter.
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018-2019 zwlab
 *
 * Version: 1.0
 *
 * Date: 12/2018
 *
 * History:
 *     12/2018: Original version created by Luqian Jiang
 *     11/2019: Luqian Jiang:
 *         change the check_assignment function for combining with modeling
 *
 ************************************************************************/
#include "pre_interface_struct.h"
#include "pre_para_assign.h"

/* We use the scan-line method to achieve assignment */
void Polygon_Fill(int IMAGE_BOT, int IMAGE_TOP, int IMAGE_LEFT, int IMAGE_RIGHT, int nx, int nz,
    int polyCorners, float *polyX, float *polyY, float value, float *array)
{
    //  public-domain code by Darel Rex Finley, 2007

    int  nodes, pixelX, pixelY, i, j;
    //int IMAGE_BOT = 0, IMAGE_TOP = nz-1, IMAGE_LEFT = 0, IMAGE_RIGHT = nx-1;
    float *nodeX = NULL, swap;

    nodeX = (float*)malloc(MAX_POLY_CORNERS*sizeof(float));

    //  Loop through the rows of the image.
    for (pixelY = IMAGE_BOT; pixelY <= IMAGE_TOP; pixelY++) {
        /* Build a list of nodes. */
        nodes = 0;
        j = polyCorners-1;
        for (i = 0; i < polyCorners; i++) {
            if (  (   polyY[i] <= (float) pixelY && polyY[j] >= (float) pixelY
                   || polyY[j] <= (float) pixelY && polyY[i] >= (float) pixelY)
                && fabs(polyY[i]-polyY[j]) > 1e-5 ) {
                nodeX[nodes++] = polyX[i]+(pixelY-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i]);
            }
            j = i;
        }

        /* Sort the nodes, via a simple "Bubble" sort. */
        i = 0;
        while (i < nodes-1) {
            if (nodeX[i] > nodeX[i+1]) {
                swap = nodeX[i];
                nodeX[i] = nodeX[i+1];
                nodeX[i+1] = swap;
                if (i)
                    i--;
            }
            else {
                i++;
            }
        }

        /* Fill the pixels between node pairs. */
        for (i = 0; i < nodes; i += 2) {
            if   (nodeX[i  ] >= IMAGE_RIGHT)
                break;
            if   (nodeX[i+1] >  IMAGE_LEFT ) {
                if (nodeX[i  ] < IMAGE_LEFT )
                    nodeX[i  ] = IMAGE_LEFT ;
                if (nodeX[i+1] > IMAGE_RIGHT)
                    nodeX[i+1] = IMAGE_RIGHT;
                for (pixelX = ceil(nodeX[i]); pixelX <= floor(nodeX[i+1]); pixelX++) {
                    array[pixelY*nx + pixelX] = value;
                }
            }
        }
    }
}


/* assignment the multi-layer model with a set of values */
int parameter_assignment(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *value, float *data,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    float *interface_x_index = NULL, *interface_z_index = NULL;
    int i_interface, i_point, ierr = 0;

    for (i_interface = 0; i_interface < number_of_interfaces; i_interface++) {

        interface_x_index = (float*) malloc(interface[i_interface].npoints_interfaces*sizeof(float));
        interface_z_index = (float*) malloc(interface[i_interface].npoints_interfaces*sizeof(float));

        /* The index of layer in coordinate */
        for (i_point = 0; i_point < interface[i_interface].npoints_interfaces; i_point++) {
            interface_x_index[i_point] = (interface[i_interface].x_loc[i_point] - xmin)/dx;
            interface_z_index[i_point] = (interface[i_interface].z_loc[i_point] - zmin)/dz;
        }

        Polygon_Fill(0, nz-1, 0, nx-1, nx, nz,
            interface[i_interface].npoints_interfaces, interface_x_index, interface_z_index,
            value[ interface[i_interface].num_of_material-1 ], data);

        free(interface_x_index);
        free(interface_z_index);
    }

    ierr = check_assignment(data, nz, nx, nghost_x1, nghost_x2, nghost_z1, nghost_z2);

    if (ierr == 0)
        return 0;
    else
        return 1;
}

/* check if all grid points are assigned */
int check_assignment(float *data, int nz, int nx,
     int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2)
{
    int iz, ix;
    for (iz = nghost_z1; iz < nz-nghost_z2; iz++) {
        for (ix = nghost_x1; ix < nx-nghost_x2; ix++) {
            if (data[nx*iz + ix] <= 0.0) {
                printf("Assignment error: row %d column %d has the wrong value, please check the interface file!\n",iz,ix);
                return 1;
            }
        }
    }
    return 0;
}
