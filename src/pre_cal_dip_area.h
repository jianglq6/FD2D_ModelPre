#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define  MAX_LAYER_POINTS 10000

/* 20190731 Try to solve the precision problem */
struct grid_id_in_layer
{
  int ifxvec;
  int ifzvec;
  float xvec;
  float zvec;
  float id;
};

void getLayer_of_grid(float *xvec, float *zvec, float *layer_x, float *layer_z,
                      int nx, int nz, int npoints_interfaces,
                      float *layer_of_grid_x, float *layer_of_grid_z, int *npoints_layer);


void cal_dip_area(float *xvec, float *zvec, float dx, float dz, int nx, int nz,
                  float *layer_x, float *layer_z, int npoints_interfaces, int *npoints_para,
                  int *ix, int *iz, int *loc_type, float *area1, float *area2, float *theta);

/* Shell sort */
void shellsort(struct grid_id_in_layer *arr, int N);

/* Merge two sorted array and Deduplication */
void merge(struct grid_id_in_layer *arr1, struct grid_id_in_layer *arr2,
           int n1, int n2, struct grid_id_in_layer *arr3, int *n3);

int findNearestNeighborIndex( float value, float *x, int len);

// Lagrange interpolation
void interp1(float *x, int x_len, float *v, float *xq, int xq_len, float *vq);
