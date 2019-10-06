#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define KAISERB 4.14
#define NSRCEXT 4      /* Kaiser window */

#define IMPULSE_DIRECT   0
#define IMPULSE_AVERAGE  1
#define IMPULSE_SINC     2


struct Src
{
    int number_of_src;
    int *stf_type_id;
    float *stf_timefactor; /* gauss t0; ricker t0; bell starting */
    float *stf_freqfactor; /* gauss a;  ricker fc; bell width    */
    float *xs;
    float *zs;
    float *Fx;
    float *Fz;
    float *Mxx;
    float *Mxz;
    float *Mzz;
};

int src_tti_moment(float *hTxx_1, float *hTxx_2,
                   float *hTxz_1, float *hTxz_2,
                   float *hTzz_1, float *hTzz_2,
                   float *xvec1, float *zvec1,
                   float *xvec2, float *zvec2,
                   int source_impulse_method, struct Src src,
                   float current_time, float dt,
                   float dx, float dz, int nx);

int src_tti_force(float *hVx_1, float *hVz_1,
                  float *hVx_2, float *hVz_2,
                  float *xvec1, float *zvec1,  /* for the integer grids  */
                  float *xvec2, float *zvec2,  /* for the half grids     */
                  int source_impulse_method, struct Src src,
                  float current_time, float dt,
                  float dx, float dz, int nx,
                  float *B01, float *B10);

float cal_source_time_function(int flag_stf_type, float t, float t0, float f0);

float calculate_windowed_impulse(float x, float b, float r);
