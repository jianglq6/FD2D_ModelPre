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


FILE *gfopen(char *filename, char *mode);

/* removes the whitespace from both ends of a string */
char *trim(char *str);

void para_read(FILE *fid,  char *string_read_value, size_t string_read_len, char *name, size_t name_len, int *ierr);

void para_read_nextparam(FILE *fid, char *string_read, size_t string_read_len, char *name, size_t name_len, int *ierr);

/* Read integer data in the Par_file */
void read_value_int(FILE *fid, char *para_name, int *value_to_read, int *ierr);

/* Read float data in the Par_file */
void read_value_float(FILE *fid, char *para_name, float *value_to_read, int *ierr);

/* Read string data in the Par_file */
void read_value_string(FILE *fid, char *para_name, char *value_to_read, int *ierr);

/* Read Boolean data in the Par_file */
void read_value_bool(FILE *fid, char *para_name, bool *value_to_read, int *ierr);

void read_value_float_next_p(FILE *fp, char *para_name, float *value_to_read, int *ierr);

void read_value_int_next_p(FILE *fp, char *para_name, int *value_to_read, int *ierr);

