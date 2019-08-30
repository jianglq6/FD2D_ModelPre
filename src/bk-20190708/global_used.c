#include "global_used.h"

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

/* Merge two sorted array and Deduplication */
void merge(float *arr1, float *arr2, int n1, int n2, float *arr3, int *n3)
{
    int i = 0, j = 0, k = 0;

    // Traverse both array
    while (i < n1 && j < n2) {

        if (arr1[i] < arr2[j]) {
            arr3[k] = arr1[i];
            k++;
            i++;
        }
        else if (arr1[i] == arr2[j]) {
            arr3[k] = arr2[j];
            k++;
            j++;
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
        if (indiceEnVector != -1)
            vq[i] = slope[indiceEnVector] * xq[i] + intercept[indiceEnVector];
        else
            vq[i] = DBL_MAX;
    }

    free(slope);
    free(intercept);

}
