/* March 27, 2020 - Andrés Nicolás Ruiz */

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
  int    i;
  double dplus,dnorm,fomega;
  double aini,aend,da,a,z;
  double omega_r,omega_m,omega_k,omega_de;
  FILE   *fd;
  	
  if (argc < 2) {
      fprintf(stdout, "\n Error. Missing input file.\n");
      fprintf(stdout, "./dplus.x <input_param>\n\n");
      exit(EXIT_FAILURE);
  }

  read_inputfile(argv[1]);

  if (UseTab_wEOS == 1) set_tables();

  aini = log10(1.0/(1.0 + InitialRedshift));
  aend = log10(1.0/(1.0 + FinalRedshift));
  da = (aend - aini)/(double)Nbins;
 
  growth(1.0, &dnorm, &fomega);	

  fd = fopen(OutputFile,"w");

  for (i=0; i<Nbins; i++) {
      	
      a = aini + da*(double)(i+1);
      a = pow(10.0,a);  
      z = 1.0/a - 1.0;
	  
      growth(a, &dplus, &fomega);	  
      omegas_a(a, &omega_r, &omega_m, &omega_k, &omega_de);

      fprintf(fd,"%f %f %f %f %f %f %f %f %f \n",a,z,dplus/dnorm,
              fomega,omega_r,omega_m,omega_k,omega_de,hubble_a(a));

  }

  fclose(fd);

}

