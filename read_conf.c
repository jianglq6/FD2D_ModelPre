/*************************************************************************
 *
 * This function is used to read the configuration file
 *
 * Authors: Luqian Jiang   Email: jlq0608@gmail.com
 *          Wei Zhang      Email: zhangwei@sustc.edu.cn
 *
 * Copyright (c) 2018-2019 zwlab
 *
 * Version: 1.0
 *
 * Date: 12/2018
 *
 * History:
 *     12/2018: Original version created by Luqian Jiang
 *
 *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "read_conf.h"


/* gfopen: elegant version of fopen() */
FILE *gfopen(char *filename, char *mode)
{
    FILE *fp;
    if ((fp = fopen(filename,mode)) == NULL) {
        fprintf(stderr, "Cannot open %s, please check your file path and run-directory. \n",filename);
        exit(1);
    }
    return fp;
}

/* removes the whitespace from both ends of a string */
char *trim(char *str)
{
  char *end;

  /* Trim leading space */
  while(isblank((unsigned char)*str)) ++str;


  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end)) end--;

  /* Write new null terminator character */
  end[1] = '\0';

  return str;
}

/*
 * Read the Config File
 *
 * Basically, the lines of the parameter file should be of the form
 * 'parameter = value', optionally followed by a #-delimited comment.
 * 'value' can be any number of space- or tab-separated words.
 * Blank lines, lines containing only white space and lines whose first
 * non-whitespace character is '#' are ignored.
 *
 * We use regular expression for parsing lines from parameter/configure file,
 * so if both parameter and value are not specified, the line is ignored.
*/
void parse_config(char *config_file, char *name, char *string_read_value, int *ierr, size_t string_read_len)
{

    FILE *file  = gfopen(config_file, "r");

    char buf[MAX_BUF_LEN];
    char key[MAX_KEY_LEN] = {0}, val[MAX_VAL_LEN] = {0};

    char pattern[] = "^[ \t]*([^# \t]+)[ \t]*=[ \t]*([^# \t]+([ \t]+[^# \t]+)*)";
    regex_t compiled_pattern; /**/
    int status;  /**/
    int regret;  /**/
    regmatch_t parameter[3];
    char *name_trim;
    char *keyword;
    char *value;
    size_t value_len;


    /* Trim the keyword name we're looking for */
    name_trim = trim(name);

    /* Compile the regular expression. */
    status = regcomp(&compiled_pattern, pattern, REG_EXTENDED);
    if (status != 0) {
        printf("regcomp return error %d\n", status);
    }


    /* Read every line in the file. */
    while(fgets(buf, MAX_BUF_LEN, file) != NULL)
    {

        /* Trim the line first */
        char *line = trim(buf);

        /* ignore and skip the comment line and blank line */
        if (line[0] == COMMENT_CHAR || line[0] == '\n')
            continue;

        /*Test if line matches the regular expression pattern,
         * if so, return position of keyword and value.
         *TBD: Maybe using strstr is much easier?           */
        regret = regexec(&compiled_pattern, line, 3, parameter, 0);

        /* If no match, check the next line. */
        if (regret == REG_NOMATCH) {
            continue;
        }

        /* If any error, bail out with an error message */
        if (regret != 0) {
            printf("regexec returned error %d\n", regret );
            *ierr = 1;
            regfree(&compiled_pattern);
            return;
        }

        /* If we have a match, extract the keyword from the line. */
        keyword = strndup(line + parameter[1].rm_so, parameter[1].rm_eo - parameter[1].rm_so);
        /* If the keyword is not the one we're looking for, check the next line. */
        if (strcmp(keyword, name_trim) != 0) {
            free(keyword);
            continue;
        }
        free(keyword);
        regfree(&compiled_pattern);
        /* If it is matches, extract the value form the line. */
        value = strndup(line + parameter[2].rm_so, parameter[2].rm_eo - parameter[2].rm_so);

        /* Clear out the return string with blanks, copy the value into it, and return. */
        memset(string_read_value, 0, string_read_len);
        value_len  = strlen(value);
        if ( value_len > string_read_len ) {
            value_len = string_read_len;
        }
        strncpy(string_read_value, value, value_len);

        trim(string_read_value);
        free(value);

        /* If it is matches, return to the superior call immediately */
        *ierr = 0;
        return ;
    }

    /* Free memory for regular expression */
    regfree(&compiled_pattern);

    /* If no keywords matches, set the error flag */
    *ierr = 1;
    fclose(file);
    return;
}


/* Read integer data in the Par_file */
void read_conf_value_int(char *config_file, char *para_name, int *value_to_read)
{
    int int_err = 0;
    int *ierr = &int_err;
    char string_read[MAX_VAL_LEN];

    parse_config(config_file, para_name, string_read, ierr ,sizeof(string_read));

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    *value_to_read = atoi(string_read);

    return;
}

/* Read float data in the Par_file */
void read_conf_value_float(char *config_file, char *para_name, float *value_to_read)
{
    int int_err = 0;
    int *ierr = &int_err;
    char string_read[MAX_VAL_LEN];

    parse_config(config_file, para_name, string_read, ierr, sizeof(string_read));

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    *value_to_read = atof(string_read);

    return;
}

/* Read string data in the Par_file */
void read_conf_value_string(char *config_file, char *para_name, char *value_to_read)
{
    int int_err = 0;
    int *ierr = &int_err;
    char string_read[MAX_VAL_LEN];

    parse_config(config_file, para_name, string_read, ierr, sizeof(string_read));

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    memset(value_to_read, 0, strlen(value_to_read));
    strncpy(value_to_read, string_read, strlen(string_read));
    trim(string_read);

    return;
}

/* Read Boolean data in the Par_file */
void read_conf_value_bool(char *config_file, char *para_name, bool *value_to_read)
{
    int int_err = 0;
    int *ierr = &int_err;
    char string_read[MAX_VAL_LEN];

    parse_config(config_file, para_name, string_read, ierr, sizeof(string_read));

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    if (strcmp(string_read, "true") == 0)
        *value_to_read = 1;
    if (strcmp(string_read, "false") == 0)
        *value_to_read = 0;

    return;
}


/* read material table, TBD: needed to be improved!!! */
/* Using structure Array?? */
void read_material(int nbmodel, char *materialfile, int *material_type, int *model_number,
    float *val1, float *val2, float *val3, float *val4, float *val5, float *val6, float *val7,
    float *val8, float *val9, float *val10, float *val11, float *val12, float *val13, int *ierr)
{
    FILE  *file  = gfopen(materialfile, "r");
    char  buf[MAX_BUF_LEN];
    int   i_split = 0;
    int   col = 15;  /* the column of material is 13, 13 + 2 */
    int   i_useful_line = 0;
    char  *split_char;
    float numtmp[col];

    /* Read every line in the file. */
    while(fgets(buf, MAX_BUF_LEN, file) != NULL)
    {
        /* Trim the line first */
        char *line = trim(buf);

        /* ignore and skip the comment line and blank line */
        if (line[0] == COMMENT_CHAR || line[0] == '\n')
            continue;

        /* Split the useful data and assignment. */
        split_char = strtok(line," ");
        i_split = 0;
        while(split_char != NULL) {
            numtmp[i_split] = atof(split_char);
            i_split ++;
            split_char = strtok(NULL," ");
        }

        /* If it has the incomplete parameters, return to the superior call */
        if (i_split != col) {
            printf("Materials file error!\n");
            *ierr = 1;
            return;
        }

        /* Assignment */
        *(model_number  + i_useful_line)  =  (int) numtmp[0];
        *(material_type + i_useful_line)  =  (int) numtmp[1];
        *(val1          + i_useful_line)  =        numtmp[2];
        *(val2          + i_useful_line)  =        numtmp[3];
        *(val3          + i_useful_line)  =        numtmp[4];
        *(val4          + i_useful_line)  =        numtmp[5];
        *(val5          + i_useful_line)  =        numtmp[6];
        *(val6          + i_useful_line)  =        numtmp[7];
        *(val7          + i_useful_line)  =        numtmp[8];
        *(val8          + i_useful_line)  =        numtmp[9];
        *(val9          + i_useful_line)  =        numtmp[10];
        *(val10         + i_useful_line)  =        numtmp[11];
        *(val11         + i_useful_line)  =        numtmp[12];
        *(val12         + i_useful_line)  =        numtmp[13];
        *(val13         + i_useful_line)  =        numtmp[14];


        i_useful_line ++;
    }

    fclose(file);
    *ierr = 0;
    return;
}

/* Read interfaces data */
/* TBD: needed to be improved!!! some error reporting should be added, not strong enough !!!*/
struct interfaces * read_interfaces(char *interfacesfile, int *number_of_interfaces)
{
    int   i_interface, npoints_interfaces, i_point;
    FILE *file = gfopen(interfacesfile, "r");
    FILE *tmp_file;
    char line[MAX_BUF_LEN];
    tmp_file = tmpfile();

    /* Read every line in the file. */
    while(fgets(line, MAX_BUF_LEN, file) != NULL)
    {
        /* ignore and skip the comment line and blank line */
        if (line[0] == COMMENT_CHAR || line[0] == '\n')
            continue;

        /* Write the useful data to the temporary file */
        fputs(line,tmp_file);
    }

    
    /* Sets the file pointer at the beginning of the stream*/
    rewind(tmp_file);

    /* Read the temporary data file, Can not use the feof() */
    while(feof(tmp_file) != EOF)
    {
        fscanf(tmp_file,"%d", number_of_interfaces);

        /* Initialize the structure array interface */
        struct interfaces *interface = (struct interfaces *)malloc(*number_of_interfaces * sizeof(struct interfaces));

        for (i_interface = 0; i_interface < *number_of_interfaces; i_interface++) {
            
            printf("%d\n",i_interface);
            fscanf(tmp_file,"%d", &npoints_interfaces);
            interface[i_interface].npoints_interfaces = npoints_interfaces;
            interface[i_interface].x_loc = (float *) malloc(npoints_interfaces*sizeof(float));
            interface[i_interface].z_loc = (float *) malloc(npoints_interfaces*sizeof(float));

            if (npoints_interfaces < 2) {
                fprintf(stderr, "Not enough interface points (minimum is 2)!\n");
                return;
            }

            /* Read interface coordinate value */
            // TBD: need to add error reporting, if it is not two column?
            for (i_point = 0; i_point < npoints_interfaces; i_point++) {
                fscanf(tmp_file,"%f", &interface[i_interface].x_loc[i_point]);
                fscanf(tmp_file,"%f", &interface[i_interface].z_loc[i_point]);
            }
            fscanf(tmp_file,"%d", &interface[i_interface].num_of_material);
        }

        /* Release the structure array interface */
        return(interface);
        free(interface);
        break;
    }
}









