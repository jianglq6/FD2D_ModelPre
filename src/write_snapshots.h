#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#define SNAP_DISPL 1
#define SNAP_VELOC 2
#define SNAP_ACCEL 3
#define SNAP_STRESS 4

int write_snapshots(float *Vx_1,  float *Vx_2,
                    float *Vz_1,  float *Vz_2,
                    float *Txx_1, float *Tzz_1,
                    float *Txz_1, float *Txx_2,
                    float *Tzz_2, float *Txz_2,
                    float *xvec1, float *xvec2,
                    float *zvec1, float *zvec2,
                    int it, int nt, float dt,   // it = it+1
                    int nx, int nz, float dx, float dz, int half_fd_stencil,
                    int *boundary_layer_number,
                    int NSTEP_BETWEEN_OUTPUT_IMAGES,
                    bool output_postscript_snapshot,
                    int imagetype_postscript, bool meshvect, bool modelvect, bool boundvect,
                    float sizemax_arrows, bool US_LETTER,
                    bool output_wavefield_dumps, int imagetype_wavefield_dumps,
                    bool use_binary_for_wavefield_dumps);

int write_wavefield_dumps(float *Vx_1,  float *Vx_2,
                          float *Vz_1,  float *Vz_2,
                          float *Txx_1, float *Tzz_1,
                          float *Txz_1, float *Txx_2,
                          float *Tzz_2, float *Txz_2,
                          float *xvec1, float *xvec2,
                          float *zvec1, float *zvec2,
                          int it, int nt, float dt,   // it = it+1
                          int nx, int nz, float dx, float dz,
                          int *boundary_layer_number,
                          int imagetype_wavefield_dumps,
                          bool use_binary_for_wavefield_dumps);

/* gfopen: elegant version of fopen() */
FILE *gfopen(char *filename, char *mode);

int write_wavefield_dumps_ssg(float *Vx, float *Vz,
                          float *Txx, float *Tzz, float *Txz,
                          float *xvec1, float *xvec2,
                          float *zvec1, float *zvec2,
                          int it, int nt, float dt,   // it = it+1
                          int nx, int nz, float dx, float dz,
                          int *boundary_layer_number,
                          int imagetype_wavefield_dumps,
                          bool use_binary_for_wavefield_dumps);

int write_snapshots_ssg(float *Vx, float *Vz,
                        float *Txx, float *Tzz, float *Txz,
                        float *xvec1, float *xvec2,
                        float *zvec1, float *zvec2,
                        int it, int nt, float dt,   // it = it+1
                        int nx, int nz, float dx, float dz, int half_fd_stencil,
                        int *boundary_layer_number,
                        int NSTEP_BETWEEN_OUTPUT_IMAGES,
                        bool output_postscript_snapshot,
                        int imagetype_postscript, bool meshvect, bool modelvect, bool boundvect,
                        float sizemax_arrows, bool US_LETTER,
                        bool output_wavefield_dumps, int imagetype_wavefield_dumps,
                        bool use_binary_for_wavefield_dumps);
