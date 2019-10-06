#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <regex.h>

#define COMMENT_CHAR '#'
#define MAX_BUF_LEN 1024
#define MAX_KEY_LEN 128
#define MAX_VAL_LEN 256

/* gfopen: elegant version of fopen() */
FILE *gfopen(char *filename, char *mode);

/* removes the whitespace from both ends of a string */
char *trim(char *str);

/* Read the Config File */
void parse_config(char *config_file, char *name, char *string_read_value, int *ierr, size_t );

/* Read integer data in the Par_file */
void read_conf_value_int(char *config_file, char *para_name, int *value_to_read);

/* Read float data in the Par_file */
void read_conf_value_float(char *config_file, char *para_name, float *value_to_read);

/* Read string data in the Par_file */
void read_conf_value_string(char *config_file, char *para_name, char *value_to_read);

/* Read Boolean data in the Par_file */
void read_conf_value_bool(char *config_file, char *para_name, bool *value_to_read);


/* Read material table */
void read_material(int nbmodel, char *materialfile, int *material_type, int *model_number,
    float *val1, float *val2, float *val3, float *val4, float *val5, float *val6, float *val7,
    float *val8, float *val9, float *val10, float *val11, float *val12, float *val13, int *ierr);

/* Read interfaces data */
struct interfaces *read_interfaces(char *interfacesfile, int *number_of_interfaces);

void read_value_float_next_p(FILE *fp, char *para_name, float *value_to_read);

void read_value_int_next_p(FILE *fp, char *para_name, int *value_to_read);

void param_read_nextparam(FILE *fid, char *string_read, size_t string_read_len, char *name, size_t name_len, int *ierr);

