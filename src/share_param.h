#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#define STATION_FILE "./DATA/STATIONS"

int get_config_info(char *config_file, int *nt, float *dt,
    int *half_fd_stencil, int *spatial_difference_method,
    float *xmin, float *dx, int *nx, float *zmin, float *dz, int *nz,
    int *filter_method,
    int *source_impulse_method, struct Src *src,
    int *seismotype, int *NSTEP_BETWEEN_OUTPUT_SEISMOS,
    bool *save_ASCII_seismograms, bool *save_binary_seismograms,
    int *nreceiver, float **xr, float **zr,
    int *boundary_type, int *boundary_layer_number,
    int *NSTEP_BETWEEN_OUTPUT_INFO,
    int *NSTEP_BETWEEN_OUTPUT_IMAGES, bool *output_postscript_snapshot,
    int *imagetype_postscript, bool *meshvect, bool *modelvect, bool *boundvect,
    float *sizemax_arrows, bool *US_LETTER,
    bool *output_wavefield_dumps, int *imagetype_wavefield_dumps,
    bool *use_binary_for_wavefield_dumps);

int write_station_coor_file(int nreceivers, float *xr, float *zr);

void getFloatMemory(float **p, int size);
