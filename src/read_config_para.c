#include "read_config_para.h"

/* gfopen: elegant version of fopen() */
FILE *gfopen(char *filename, char *mode)
{
    FILE *fp;
    if ((fp = fopen(filename,mode)) == NULL) {
        fprintf(stderr, "Cannot open %s\n, please check your file path and run-directory. \n",filename);
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
void para_read(FILE *fid, char *string_read_value, size_t string_read_len, char *name, size_t name_len, int *ierr)
{

    char buf[MAX_BUF_LEN];
    char key[MAX_KEY_LEN] = {0}, val[MAX_VAL_LEN] = {0};

    char pattern[] = "^[ \t]*([^# \t]+)[ \t]*=[ \t]*([^# \t]+([ \t]+[^# \t]+)*)";
    regex_t compiled_pattern; /**/
    int status;  /**/
    int regret;  /**/
    regmatch_t parameter[3];
//    char *name_trim;
    char *keyword;
    char *value;
    size_t value_len;


    /* Trim the keyword name we're looking for */
//    name_trim = trim(name);

    /* Compile the regular expression. */
    status = regcomp(&compiled_pattern, pattern, REG_EXTENDED);
    if (status != 0) {
        printf("regcomp return error %d\n", status);
    }

    /* Position the open file to the begining */
    if (fseek(fid, 0, SEEK_SET) !=0 ) {
        printf("Can't seek to beginning of parameter file \n");
        *ierr = 1;
        regfree(&compiled_pattern);
        return;
    }


    /* Read every line in the file. */
    while(fgets(buf, MAX_BUF_LEN, fid) != NULL)
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
  //      if (strcmp(keyword, name_trim) != 0) {
        if (strcmp(keyword, name) != 0) {
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
//    fclose(file);
    return;
}


void para_read_nextparam(FILE *fid, char *string_read_value, size_t string_read_len, char *name, size_t name_len, int *ierr)
{

    char buf[MAX_BUF_LEN];
    char key[MAX_KEY_LEN] = {0}, val[MAX_VAL_LEN] = {0};

    char pattern[] = "^[ \t]*([^# \t]+)[ \t]*=[ \t]*([^# \t]+([ \t]+[^# \t]+)*)";
    regex_t compiled_pattern; /**/
    int status;  /**/
    int regret;  /**/
    regmatch_t parameter[3];
//    char *name_trim;
    char *keyword;
    char *value;
    size_t value_len;


    /* Trim the keyword name we're looking for */
//    name_trim = trim(name);

    /* Compile the regular expression. */
    status = regcomp(&compiled_pattern, pattern, REG_EXTENDED);
    if (status != 0) {
        printf("regcomp return error %d\n", status);
    }

    // Do not reposition to BOF

    /* Read next lines in the file until we have a match */
    while(fgets(buf, MAX_BUF_LEN, fid) != NULL)
    {

        /* Trim the line first */
        char *line = trim(buf);

        /* ignore and skip the comment line and blank line */
        if (line[0] == COMMENT_CHAR || line[0] == '\n')
            continue;
        //if (line[strlen(line)-1] = '\n')
        //    line[strlen(line)-1] = '\0';

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
  //      if (strcmp(keyword, name_trim) != 0) {
        if (strcmp(keyword, name) != 0) {
            printf("keyword returned wrong parameter %s instead of %s \n", keyword,name);
            free(keyword);
            *ierr = 1;
            regfree(&compiled_pattern);
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
//    fclose(file);
    return;
}

/* Read integer data in the Par_file */
void read_value_int(FILE *fid, char *para_name, int *value_to_read, int *ierr)
{
    char string_read[MAX_VAL_LEN];

    para_read(fid,  string_read, sizeof(string_read), para_name, sizeof(para_name), ierr);

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    *value_to_read = atoi(string_read);

    return;
}

/* Read float data in the Par_file */
void read_value_float(FILE *fid, char *para_name, float *value_to_read, int *ierr)
{

    char string_read[MAX_VAL_LEN];

    para_read(fid, string_read, sizeof(string_read), para_name, sizeof(para_name), ierr);

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    *value_to_read = atof(string_read);

    return;
}

/* Read string data in the Par_file */
void read_value_string(FILE *fid, char *para_name, char *value_to_read, int *ierr)
{

    char string_read[MAX_VAL_LEN];

    para_read(fid, string_read, sizeof(string_read), para_name, sizeof(para_name), ierr);

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    memset(value_to_read, 0, strlen(value_to_read));
    strncpy(value_to_read, string_read, strlen(string_read));
    trim(string_read);

    return;
}

/* Read Boolean data in the Par_file */
void read_value_bool(FILE *fid, char *para_name, bool *value_to_read, int *ierr)
{

    char string_read[MAX_VAL_LEN];

    para_read(fid, string_read, sizeof(string_read), para_name, sizeof(para_name), ierr);

    /* If read nothing, return to the superior call, then use the default parameter. */
    if (*ierr != 0) return;

    if (strcmp(string_read, "true") == 0)
        *value_to_read = 1;
    if (strcmp(string_read, "false") == 0)
        *value_to_read = 0;

    return;
}

void read_value_float_next_p(FILE *fp, char *para_name, float *value_to_read, int *ierr)
{
    char string_read[MAX_VAL_LEN];

    para_read_nextparam(fp, string_read, sizeof(string_read), para_name, sizeof(para_name), ierr);

    if (*ierr != 0) return;

    *value_to_read = atof(string_read);

    return;
}

void read_value_int_next_p(FILE *fp, char *para_name, int *value_to_read, int *ierr)
{
    char string_read[MAX_VAL_LEN];

    para_read_nextparam(fp, string_read, sizeof(string_read), para_name, sizeof(para_name), ierr);

    if (*ierr != 0) return;

    *value_to_read = atoi(string_read);

    return;
}
