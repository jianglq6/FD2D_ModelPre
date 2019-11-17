#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define KAISERB 4.14
#define NSRCEXT 4      /* Kaiser window */

#define IMPULSE_DIRECT   1
#define IMPULSE_AVERAGE  2
#define IMPULSE_SINC     3


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

int src_moment_lebedev(float *hTxx_1, float *hTxx_2,
                       float *hTxz_1, float *hTxz_2,
                       float *hTzz_1, float *hTzz_2,
                       float *xvec1, float *zvec1,
                       float *xvec2, float *zvec2,
                       int source_impulse_method, struct Src src,
                       float current_time, float dt,
                       float dx, float dz, int nx);

int src_force_lebedev(float *hVx_1, float *hVz_1,
                      float *hVx_2, float *hVz_2,
                      float *xvec1, float *zvec1,  /* for the integer grids  */
                      float *xvec2, float *zvec2,  /* for the half grids     */
                      int source_impulse_method, struct Src src,
                      float current_time, float dt,
                      float dx, float dz, int nx);

/* Moment source (sinc distribution) */
int src_sinc_moment_lebedev(float *hTxx_1, float *hTxx_2,
                            float *hTxz_1, float *hTxz_2,
                            float *hTzz_1, float *hTzz_2,
                            float *xvec1, float *zvec1,
                            float *xvec2, float *zvec2,
                            struct Src src,
                            float current_time, float dt,
                            float dx, float dz, int nx);

/* Force source (sinc distribution) */
int src_sinc_force_lebedev(float *hVx_1, float *hVz_1,
                           float *hVx_2, float *hVz_2,
                           float *xvec1, float *zvec1,  /* for the integer grids  */
                           float *xvec2, float *zvec2,  /* for the half grids     */
                           struct Src src,
                           float current_time, float dt,
                           float dx, float dz, int nx);

float cal_source_time_function(int flag_stf_type, float t, float t0, float f0);

double calculate_windowed_impulse(float x, float b, float r);

/* for staggered scheme */
int src_moment_ssg(float *hTxx, float *hTxz, float *hTzz,
                         float *xvec1, float *zvec1,
                         float *xvec2, float *zvec2,
                         int source_impulse_method, struct Src src,
                         float current_time, float dt,
                         float dx, float dz, int nx);
/* Force source */
int src_force_ssg(float *hVx, float *hVz,
                        float *xvec1, float *zvec1,  /* for the integer grids  */
                        float *xvec2, float *zvec2,  /* for the half grids     */
                        int source_impulse_method, struct Src src,
                        float current_time, float dt,
                        float dx, float dz, int nx);

/* Moment source */
int src_sinc_moment_ssg(float *hTxx, float *hTxz, float *hTzz,
                              float *xvec1, float *zvec1,
                              float *xvec2, float *zvec2,
                              struct Src src,
                              float current_time, float dt,
                              float dx, float dz, int nx);

/* Force source */
int src_sinc_force_ssg(float *hVx, float *hVz,
                             float *xvec1, float *zvec1,  /* for the integer grids  */
                             float *xvec2, float *zvec2,  /* for the half grids     */
                             struct Src src,
                             float current_time, float dt,
                             float dx, float dz, int nx);
