#include "elastic2d_src.h"
#include "read_config_para.h"
#include "share_param.h"

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
    bool *use_binary_for_wavefield_dumps)

{

    int is, ireceiverset, nreceiversets, i, ireceiver=0, nreceivers = 0;
    int *nrec=NULL, ierr = 0;
    float *xdeb=NULL, *zdeb=NULL, *xfin=NULL, *zfin=NULL, *tmp=NULL;
    bool use_existing_station;
    FILE *fid=gfopen(config_file, "r");
    /* For matlab plot */
    FILE *fp_mfile = gfopen("./mfiles/configure","w");

    fprintf(stdout, "===========================================\n");
    fprintf(stdout, "       simulation input parameters         \n");
    fprintf(stdout, "===========================================\n");

    /*========================= time information ==========================*/
    read_value_int(fid, "NSTEP", nt,&ierr);
    read_value_float(fid, "DT",  dt, &ierr);
    read_value_int(fid,"half_fd_stencil", half_fd_stencil, &ierr);
    read_value_int(fid,"spatial_difference_method", spatial_difference_method, &ierr);
    read_value_int(fid,"filter_method", filter_method, &ierr);
    fprintf(stdout, "\n-------\nTime information \n");
    fprintf(stdout, "   Time step DT = %f, total number of time steps NSTEP = %d.\n",
             *dt, *nt);
    fprintf(fp_mfile, "%-42.6f    dt\n", *dt);
    fprintf(fp_mfile, "%-42d    nt\n", *nt);

    /*========================= grid information ==========================*/
    read_value_float(fid, "xmin", xmin, &ierr);
    read_value_float(fid, "zmin", zmin, &ierr);
    read_value_float(fid, "dx"  , dx, &ierr);
    read_value_float(fid, "dz"  , dz, &ierr);
    read_value_int(fid,   "nx"  , nx, &ierr);
    read_value_int(fid,   "nz"  , nz, &ierr);
    fprintf(stdout, "\n-------\nGeometry and the grid setting of the model.\n");
    fprintf(stdout, "   xmin = %f, dx = %f, nx = %d\n", *xmin, *dx, *nx);
    fprintf(stdout, "   zmin = %f, dz = %f, nz = %d\n", *zmin, *dz, *nz);

    fprintf(fp_mfile, "%-14.6f%-14.6f%-14d    xmin dx nx\n", *xmin, *dx, *nx);
    fprintf(fp_mfile, "%-14.6f%-14.6f%-14d    zmin dz nz\n", *zmin, *dz, *nz);

    /*========================== source information ========================*/
    //-- initiate source
    //!!TODO: Pass the structure pointer, then alloc the size of these elements??
    ///  There might be a segment fault problem??
    ///  should molloc in main(), then realloc?
    read_value_int(fid,"source_impulse_method", source_impulse_method, &ierr);
    read_value_int(fid, "NSOURCES", &(src->number_of_src), &ierr);
    src->stf_type_id    = (int*)   malloc(src->number_of_src*sizeof(int));
    src->stf_timefactor = (float*) malloc(src->number_of_src*sizeof(float));
    src->stf_freqfactor = (float*) malloc(src->number_of_src*sizeof(float));
    src->xs  = (float*) malloc(src->number_of_src*sizeof(float));
    src->zs  = (float*) malloc(src->number_of_src*sizeof(float));
    src->Fx  = (float*) malloc(src->number_of_src*sizeof(float));
    src->Fz  = (float*) malloc(src->number_of_src*sizeof(float));
    src->Mxx = (float*) malloc(src->number_of_src*sizeof(float));
    src->Mxz = (float*) malloc(src->number_of_src*sizeof(float));
    src->Mzz = (float*) malloc(src->number_of_src*sizeof(float));

    fprintf(stdout, "\n-------\nSource parameters\n");
    fprintf(stdout, "   Source impulse method: %d\n", *source_impulse_method);
    fprintf(stdout, "   number of sources: %d\n", src->number_of_src);


    for (is = 0; is < src->number_of_src; is++) {
        read_value_float_next_p(fid, "xs" ,                &(src->xs[is]            ), &ierr);
        read_value_float_next_p(fid, "zs" ,                &(src->zs[is]            ), &ierr);
        read_value_int_next_p(  fid, "time_function_type", &(src->stf_type_id[is]   ), &ierr);
        read_value_float_next_p(fid, "freqfactor" ,        &(src->stf_freqfactor[is]), &ierr);
        read_value_float_next_p(fid, "timefactor" ,        &(src->stf_timefactor[is]), &ierr);
        read_value_float_next_p(fid, "Fx" ,                &(src->Fx[is]            ), &ierr);
        read_value_float_next_p(fid, "Fz" ,                &(src->Fz[is]            ), &ierr);
        read_value_float_next_p(fid, "Mxx",                &(src->Mxx[is]           ), &ierr);
        read_value_float_next_p(fid, "Mzz",                &(src->Mzz[is]           ), &ierr);
        read_value_float_next_p(fid, "Mxz",                &(src->Mxz[is]           ), &ierr);

        fprintf(stdout,"   --- Source %d --- \n", is+1);
        fprintf(stdout,"    Source location: xs = %f zs = %f\n", src->xs[is], src->zs[is]);
        fprintf(stdout,"    source time function (second derivative of a Gaussian (a.k.a. Ricker) = 1, first derivative of a Gaussian = 2, Gaussian = 3): %d\n", src->stf_type_id[is]);
        fprintf(stdout,"    domainant source frequency (Hz): %f; time delay t0: %f\n",
            src->stf_timefactor[is], src->stf_freqfactor[is]);
        fprintf(stdout,"    Fx = %f, Fz = %f\n", src->Fx[is], src->Fz[is]);
        fprintf(stdout,"    Mxx = %f, Mzz = %f, Mxz = %f\n", src->Mxx[is], src->Mzz[is], src->Mxz[is]);
    }

    /*===================== receiver information ==========================*/
    read_value_int( fid, "seismotype", seismotype, &ierr);
    read_value_int( fid, "NSTEP_BETWEEN_OUTPUT_SEISMOS", NSTEP_BETWEEN_OUTPUT_SEISMOS, &ierr);
    read_value_bool(fid, "save_binary_seismograms", save_binary_seismograms, &ierr);
    read_value_bool(fid, "save_ASCII_seismograms", save_ASCII_seismograms, &ierr);
    read_value_bool(fid, "use_existing_STATIONS" , &use_existing_station, &ierr);
    read_value_int( fid, "nreceiversets", &nreceiversets, &ierr);
    fprintf(stdout, "\n-------\nStations\n");

    nrec = (int*)   malloc(nreceiversets*sizeof(int));
    xdeb = (float*) malloc(nreceiversets*sizeof(float));
    zdeb = (float*) malloc(nreceiversets*sizeof(float));
    xfin = (float*) malloc(nreceiversets*sizeof(float));
    zfin = (float*) malloc(nreceiversets*sizeof(float));

    if(!use_existing_station) {
        for ( ireceiverset = 0; ireceiverset < nreceiversets; ireceiverset++ ) {
            read_value_int_next_p(  fid, "nrec", nrec+ireceiverset, &ierr);
            read_value_float_next_p(fid, "xdeb", xdeb+ireceiverset, &ierr);
            read_value_float_next_p(fid, "zdeb", zdeb+ireceiverset, &ierr);
            read_value_float_next_p(fid, "xfin", xfin+ireceiverset, &ierr);
            read_value_float_next_p(fid, "zfin", zfin+ireceiverset, &ierr);
        }
    /* transfer to the station location and write to the OUTPUT/STATION file */
        for (ireceiverset = 0; ireceiverset < nreceiversets; ireceiverset++) {
            nreceivers += nrec[ireceiverset];
        }

        printf("nreceiver: %d\n",nreceivers);

        getFloatMemory(xr, nreceivers);
        getFloatMemory(zr, nreceivers);

        printf("%f\n\n",(*xr)[5]);


        for(ireceiverset = 0; ireceiverset < nreceiversets; ireceiverset++) {
            for (i = 0; i < nrec[ireceiverset]; i++) {
                (*xr)[ireceiver] = xdeb[ireceiverset] + (xfin[ireceiverset] - xdeb[ireceiverset])/
                                   (nrec[ireceiverset]-1) * i;
                (*zr)[ireceiver] = zdeb[ireceiverset] + (zfin[ireceiverset] - zdeb[ireceiverset])/
                                   (nrec[ireceiverset]-1) * i;
                ireceiver++;
            }
        }

        /* write station coordinates to DATA/STATION */
        write_station_coor_file(nreceivers, *xr, *zr);
        fprintf(stdout, "writing the %s file\n", STATION_FILE );
    }
    *nreceiver = nreceivers;

    fprintf(stdout, "   There are %6d receivers\n", nreceivers);
    fprintf(stdout, "   Target positions (x,z) of the %6d receivers\n", nreceivers);

    for (ireceiver = 0; ireceiver<nreceivers; ireceiver++)
    {
        fprintf(stdout, "   Receiver %6d = %20.6f %20.6f\n", ireceiver+1, (*xr)[ireceiver], (*zr)[ireceiver]);
    }

    free(nrec);
    free(xdeb); free(zdeb);
    free(xfin); free(zfin);

    /*======================== boundary information ==========================*/
    read_value_int(fid, "boundary_type_left",  boundary_type+0, &ierr);
    read_value_int(fid, "boundary_type_right", boundary_type+1, &ierr);
    read_value_int(fid, "boundary_type_top",   boundary_type+2, &ierr);
    read_value_int(fid, "boundary_type_bottom",boundary_type+3, &ierr);

    read_value_int(fid, "boundary_layer_number_left",  boundary_layer_number+0, &ierr);
    read_value_int(fid, "boundary_layer_number_right", boundary_layer_number+1, &ierr);
    read_value_int(fid, "boundary_layer_number_top",   boundary_layer_number+2, &ierr);
    read_value_int(fid, "boundary_layer_number_bottom",boundary_layer_number+3, &ierr);
    fprintf(stdout, "\n-------\nBoundary conditions\n");
    fprintf(stdout, "   Type of boundary condition (1. exp, 2. ade-cfs-pml 0. free-surface B.C)\n");
    fprintf(stdout, "    top:%d, bottom: %d, left: %d, right: %d\n",
            boundary_type[2], boundary_type[3], boundary_type[0], boundary_type[1]);
    fprintf(stdout, "   Number of boundary layer (if boundary type is free-surface, its number of boundary layer is ignored)\n");
    fprintf(stdout, "    top:%d, bottom: %d, left: %d, right: %d\n",
            boundary_layer_number[2], boundary_layer_number[3], boundary_layer_number[0], boundary_layer_number[1]);

    /*======================= display parameters =============================*/
    read_value_int(fid, "NSTEP_BETWEEN_OUTPUT_INFO", NSTEP_BETWEEN_OUTPUT_INFO, &ierr);


    /*====================== movies/images/snapshots =========================*/
    read_value_int(fid, "NSTEP_BETWEEN_OUTPUT_IMAGES", NSTEP_BETWEEN_OUTPUT_IMAGES, &ierr);
    read_value_int(fid, "imagetype_postscript", imagetype_postscript, &ierr);
    read_value_float(fid, "sizemax_arrows", sizemax_arrows, &ierr);
    read_value_bool(fid,"output_postscript_snapshot",  output_postscript_snapshot , &ierr);
    read_value_bool(fid,"meshvect",  meshvect, &ierr);
    read_value_bool(fid,"modelvect", modelvect, &ierr);
    read_value_bool(fid,"boundvect", boundvect, &ierr);
    read_value_bool(fid,"US_LETTER", US_LETTER, &ierr);

    read_value_bool(fid,"output_wavefield_dumps", output_wavefield_dumps, &ierr);
    read_value_bool(fid,"use_binary_for_wavefield_dumps", use_binary_for_wavefield_dumps, &ierr);
    read_value_int(fid, "imagetype_wavefield_dumps", imagetype_wavefield_dumps, &ierr);

    fclose(fid);
    fclose(fp_mfile);
    fprintf(stdout, "\n");
    return ierr;
}


/* write station coordinates to DATA/SOURCE */
int write_station_coor_file(int nreceivers, float *xr, float *zr)
{
    int ir;
    FILE *fp = NULL;
    fp = fopen(STATION_FILE, "w");
    for (ir = 0; ir < nreceivers; ir++)
        fprintf(fp, "S%04d     AA %20.6f %20.6f\n", (ir+1), xr[ir], zr[ir]);
    fclose(fp);
    return;
}

void getFloatMemory(float **p, int size)
{
    *p = malloc(sizeof(float)*size);
    if (NULL == *p) {
        printf("malloc failed! \n");
        *p = NULL;
    }
}
