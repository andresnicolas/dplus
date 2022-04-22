
#include "allvars.h"
#include "proto.h"

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

  strcpy(tag[nt],"OmegaMatter0");
  addr[nt] = &OmegaMatter0;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"OmegaRadiation0");
  addr[nt] = &OmegaRadiation0;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"OmegaCurvature0");
  addr[nt] = &OmegaCurvature0;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"OmegaDarkEnergy0");
  addr[nt] = &OmegaDarkEnergy0;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"w0EoS");
  addr[nt] = &w0EoS;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt],"waEoS");
  addr[nt] = &waEoS;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"Hubble0");
  addr[nt] = &Hubble0;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"UseTab_wEoS");
  addr[nt] = &UseTab_wEoS;
  id[nt++] = INT;

  strcpy(tag[nt],"SecondOrderDplus");
  addr[nt] = &SecondOrderDplus;
  id[nt++] = INT;

  strcpy(tag[nt],"wEoSFile");
  addr[nt] = wEoSFile;
  id[nt++] = STRING;

  strcpy(tag[nt],"OutputFile");
  addr[nt] = OutputFile;
  id[nt++] = STRING;

  strcpy(tag[nt],"InitialRedshift");
  addr[nt] = &InitialRedshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"FinalRedshift");
  addr[nt] = &FinalRedshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt],"Nbins");
  addr[nt] = &Nbins;
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

