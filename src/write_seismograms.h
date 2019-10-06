#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>


#define KAISERB 4.14
#define NRECSAMPLE 4      /* Kaiser window */


int write_seismograms(char *config_file,
                      float *Vx_1, float *Vx_2,
                      float *Vz_1, float *Vz_2,
                      float *Txx_1, float *Tzz_1,
                      float *Txz_1, float *Txx_2,
                      float *Tzz_2, float *Txz_2,
                      float *xvec1, float *xvec2,
                      float *zvec1, float *zvec2,
                      int it, int nt, float t0, float dt,
                      int nx, int nz, float dx, float dz);

int write_velocity_seismo(float *Vx_1, float *Vx_2,
                          float *Vz_1, float *Vz_2, int nx,
                          int nreceiver, float *xr, float *zr,
                          int it, int nt, float dt, float t0,
                          int NSTEP_BETWEEN_OUTPUT_SEISMOS,
                          bool save_binary_seismograms,
                          bool save_ASCII_seismograms,
                          float *xvec1, float *xvec2,
                          float *zvec1, float *zvec2,
                          float dx, float dz);

float calculate_windowed_interp(float x, float b, float r);

int read_station_info(char *config_file, int *nreceivers, float *xr, float *zr);

int write_station_coor_file(int nreceiver, float *xr, float *zr);
