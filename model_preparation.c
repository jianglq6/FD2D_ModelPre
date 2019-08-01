#include <stdio.h>
#include <stdlib.h>
#include "read_conf.h"

int main(int agrc, char *argv)
{
    
    char  config_file[] =  "../DATA/Par_file";

    /* Grid information */
    float xmin, xmax, dx, zmin, zmax, dz;
    char  xmin_key[] = "xmin", xmax_key[] = "xmax", dx_key[] = "dx", dz_key[] = "dz", zmin_key[] = "zmin", zmax_key[] = "zmax";

    int   nbmodels, effective_para_method;
    char  nbmodels_key[] = "nbmodels", effective_para_method_key[] = "effective_para_method";
    bool  use_existing_model;
    char  use_existing_model_key[] = "use_existing_model";
    char  interfacesfile[MAX_VAL_LEN], interfacesfile_key[] = "interfacesfile";
    char  materialfile[MAX_VAL_LEN], materialfile_key[] = "materialfile";

    /* Read the parameters I need */
    read_conf_value_float(config_file, xmin_key, &xmin);
    read_conf_value_float(config_file, xmax_key, &xmax);
    read_conf_value_float(config_file, dx_key  , &dx  );
    read_conf_value_float(config_file, dz_key  , &dz  );
    read_conf_value_float(config_file, zmax_key, &zmax);
    read_conf_value_float(config_file, zmin_key, &zmin);
    read_conf_value_int(config_file, nbmodels_key, &nbmodels);
    read_conf_value_int(config_file, effective_para_method_key, &effective_para_method);
    read_conf_value_bool(config_file, use_existing_model_key, &use_existing_model);

    if (!use_existing_model) {
        read_conf_value_string(config_file,  materialfile_key, materialfile);
        read_conf_value_string(config_file,interfacesfile_key,interfacesfile);
    }
        


        

}

