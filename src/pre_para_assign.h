#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define MAX_POLY_CORNERS 1000

void Polygon_Fill(int IMAGE_BOT, int IMAGE_TOP, int IMAGE_LEFT, int IMAGE_RIGHT, int nx, int nz,
    int polyCorners,float *polyX,float *polyY,float value,float *array);

int parameter_assignment(struct interfaces *interface, int number_of_interfaces,
        float xmin, float zmin, float dx, float dz, int nx, int nz, float *value, float *data,
        int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2);

/* check if all grid points are assigned */
int check_assignment(float *data, int nz, int nx,
     int nghost_x1, int nghost_x2, int nghost_z1, int nghost_z2);
