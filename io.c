
#include "io.h"

struct global_parameters P;

void read_inputfile(char *filename)
{
#define FLOAT   1
#define STRING  2
#define INT     3
#define MAXTAGS 300

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  FILE *fd;
  char buf[200],buf1[200];
  char buf2[200],buf3[200];
  int  errorFlag = 0;

  nt = 0;

  strcpy(tag[nt],"OmegaMatter");
  addr[nt] = &P.OmegaMatter;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"OmegaRadiation");
  addr[nt] = &P.OmegaRadiation;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"OmegaCurvature");
  addr[nt] = &P.OmegaCurvature;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"OmegaDarkEnergy");
  addr[nt] = &P.OmegaDarkEnergy;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"w0EoS");
  addr[nt] = &P.w0EoS;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt],"waEoS");
  addr[nt] = &P.waEoS;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"Hubble");
  addr[nt] = &P.Hubble;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"UseTab_wEoS");
  addr[nt] = &P.UseTab_wEoS;
  id[nt++] = INT;

  strcpy(tag[nt],"wEoSFile");
  addr[nt] = P.wEoSFile;
  id[nt++] = STRING;

  strcpy(tag[nt],"OutputFile");
  addr[nt] = P.OutputFile;
  id[nt++] = STRING;

  strcpy(tag[nt],"InitialRedshift");
  addr[nt] = &P.InitialRedshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"FinalRedshift");
  addr[nt] = &P.FinalRedshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"Nbins");
  addr[nt] = &P.Nbins;
  id[nt++] = INT;

  if ((fd = fopen(filename, "r"))) {
      while (!feof(fd)) {

	  buf[0] = 0;
	  fgets(buf, 200, fd);

	  if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
	  if (buf1[0] == '%') continue;

	  for (i = 0, j = -1; i < nt; i++) {
	      if (strcmp(buf1, tag[i]) == 0) {
		 j = i;
		 tag[i][0] = 0;
		 break;
	      }
	  }

	  if(j >= 0) {

	      switch (id[j]) {

		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;

		case STRING:
		  strcpy(addr[j], buf2);
		  break;

		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;

		}
	    } else {
		fprintf(stdout, "Error in file %s:  Tag '%s' not allowed or multiple defined. \n", filename, buf1);
	        errorFlag = 1;
	    }
	}
        fclose(fd);

  } else {

	fprintf(stdout, "Error. Parameter file %s not found.\n", filename);
        errorFlag = 1;

  }

  for(i = 0; i < nt; i++) {
      if(*tag[i]) {
	  fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], filename);
	  errorFlag = 1;
	}
  }

  if(errorFlag) exit(EXIT_FAILURE);

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}

