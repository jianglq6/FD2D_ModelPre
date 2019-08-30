#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* structure array of the interfaces, maybe not here */
struct interfaces{
    int   num_of_material;
    int   npoints_interfaces;
    float *x_loc;
    float *z_loc;
};

/* Shell sort */
void shellsort(float *arr, int N);

/* Merge two sorted array and Deduplication */
void merge(float *arr1, float *arr2, int n1, int n2, float *arr3, int *n3);

int findNearestNeighborIndex( float value, float *x, int len);

// Lagrange interpolation
void interp1(float *x, int x_len, float *v, float *xq, int xq_len, float *vq);
