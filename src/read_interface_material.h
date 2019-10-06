
struct interfaces * read_interfaces(char *interfacesfile, int *number_of_interfaces);

void read_material(int nbmodel, char *materialfile, int *material_type, int *model_number,
    float *val1, float *val2, float *val3, float *val4, float *val5, float *val6, float *val7,
    float *val8, float *val9, float *val10, float *val11, float *val12, float *val13, int *ierr);
