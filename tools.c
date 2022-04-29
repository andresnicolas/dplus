
#include "tools.h"

/* Count lines from an ASCII file */

int count_lines(char *filename) 
{
    FILE *fp; 
    int  count = 0;
    char string[256];
    char error_text[500];

    fp = fopen(filename, "r");
    if (fp == NULL) {
       sprintf(error_text,"File %s not found \n",filename);	   
       error(error_text); 
    }
    
    while (fgets(string,256,fp)) count++;
    fclose(fp);

    return count;
}

/* Error message and exit code */

void error(char text[])
{
  fprintf(stderr,"\nERROR!\n");
  fprintf(stderr,"%s\n",text);
  fprintf(stderr,"Aborting program.\n\n");
  exit(EXIT_FAILURE);
}

