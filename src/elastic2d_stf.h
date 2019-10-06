
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SEIS_ZERO 1e-10

//#define SIG_SVF_BELL     100
//#define SIG_SVF_TRIANGLE 110
//#define SIG_SVF_GAUSS    120
#define SIG_STF_RICKER      1
#define SIG_STF_GAUSS_DERIV 2
#define SIG_STF_GAUSS       3
//#define SIG_STF_BELL     220
//#define SIG_STF_DELTA    230
//
#define STF_FILE "./OUTPUT/source_time_function.txt"

int prepare_source_time_function(int nt, float dt, float t0, float fc, int stf_type_id, int is);

float fun_gauss(float t, float a, float t0);

float fun_gauss_deriv(float t, float a, float t0);

float fun_ricker(float t, float fc, float t0);

float fun_ricker_deriv(float t, float fc, float t0);

float fun_bshift(float t, float riset, float t0);

float fun_delta(float t, float t0);

float fun_bell(float t, float riset);

float fun_bell_deriv(float t, float riset);

float fun_bell_int(float t, float riset);

float fun_triangle(float t, float riset);

float fun_triangle_int(float t, float riset);
