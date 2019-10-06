
#include "pre_interface_struct.h"
#include "read_config_para.h"

/* Read interfaces data */
/* TBD: needed to be improved!!! some error reporting should be added, not strong enough !!!*/
struct interfaces *read_interfaces(char *interfacesfile, int *number_of_interfaces)
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

        fscanf(tmp_file, "%d", number_of_interfaces);


        /* Initialize the structure array interface */
        struct interfaces *interface = (struct interfaces *)malloc(*number_of_interfaces * sizeof(struct interfaces));

        for (i_interface = 0; i_interface < *number_of_interfaces; i_interface++) {

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

        int model_num;
        model_num                     =  (int) numtmp[0];
        *(material_type + (model_num-1))  =  (int) numtmp[1];
        *(val1          + (model_num-1))  =        numtmp[2];
        *(val2          + (model_num-1))  =        numtmp[3];
        *(val3          + (model_num-1))  =        numtmp[4];
        *(val4          + (model_num-1))  =        numtmp[5];
        *(val5          + (model_num-1))  =        numtmp[6];
        *(val6          + (model_num-1))  =        numtmp[7];
        *(val7          + (model_num-1))  =        numtmp[8];
        *(val8          + (model_num-1))  =        numtmp[9];
        *(val9          + (model_num-1))  =        numtmp[10];
        *(val10         + (model_num-1))  =        numtmp[11];
        *(val11         + (model_num-1))  =        numtmp[12];
        *(val12         + (model_num-1))  =        numtmp[13];
        *(val13         + (model_num-1))  =        numtmp[14];

        // Need a error reporting to check if model_num = 0:nbmodel-1;

        i_useful_line ++;
    }

    /* If model number */
    if (i_useful_line != nbmodel) {
        printf(" Materials file error!\n");
        *ierr = 1;
        return;
    }

    fclose(file);
    *ierr = 0;
    return;
}
