#include "utility_vector.h"

float *creat_float_array(int n, float v)
{
    int i;
    float *array;

    if (n <= 0 ) {
        fprintf(stderr, "Error: the length of array is wrong!\n");
        return NULL;
    }

    array = (float*) malloc( n * sizeof(float));

    if ( array == NULL ) {
        fprintf(stderr, "Error: while allocating!\n");
    }

    for (i = 0; i < n; i++ ) {
        array[i] = v;
    }
    return array;
}

int free_float_array(float *p)
{
   //if (p==NULL) {
   //    fprintf(stderr, "Pointer has been freed!\n");
   //    return -1;
   //}

    free(p);
    p = NULL;

    return 0;
}
