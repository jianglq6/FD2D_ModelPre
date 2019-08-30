#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define  MAX_LAYER_POINTS 10000
#define  SYS_ERR 0.0005

void getLayer_of_grid(float *xvec, float *zvec, float *layer_x, float *layer_z,
                      int nx, int nz, int npoints_interfaces,
                      float *layer_of_grid_x, float *layer_of_grid_z, int *npoints_layer);

void cal_dip_area(float *xvec, float *zvec, float dx, float dz, int nx, int nz,
                  float *layer_x, float *layer_z, int npoints_interfaces, int *npoints_para,
                  int *ix, int *iz, int *loc_type, float *area1, float *area2, float *theta);

