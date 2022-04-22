/* March 27, 2020 - Andrés Nicolás Ruiz */

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
  int    i;
  double D1,D2,f1,f2;
  double dplus,fomega,dplus_2,fomega_2;
  double aini,aend,da,a,z,dnorm;
  double omega_r,omega_m,omega_k,omega_de;
  FILE   *fd;
  	
  if (argc < 2) {
      fprintf(stdout, "\n Error. Missing input file.\n");
      fprintf(stdout, "./dplus.x <input_param>\n\n");
      exit(EXIT_FAILURE);
  }

  read_inputfile(argv[1]);

  if (UseTab_wEoS == 1) set_dark_energy_tables();
  if (SecondOrderDplus == 1) set_dplus_spline();

  growth_factor(1.0, &dnorm, &fomega);	  
  aini = log10(1.0/(1.0 + InitialRedshift));
  aend = log10(1.0/(1.0 + FinalRedshift));
  da = (aend - aini)/(double)(Nbins-1);

  fd = fopen(OutputFile,"w");
  
  for (i=0; i<Nbins; i++) {

      a = pow(10.0,aini + da*(double)i);
      z = 1.0/a - 1.0;

      density_parameters(a, &omega_r,&omega_m, &omega_k, &omega_de);
      growth_factor(a, &dplus, &fomega);
      fitting_formulas(a, &D1, &D2, &f1, &f2); 

      if (SecondOrderDplus == 1) {
	 growth_factor_2(a, &dplus_2, &fomega_2);	
         fprintf(fd,"%10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e \n",
         	      a,z,omega_r,omega_m,omega_k,omega_de,hubble_parameter(a),dplus,fomega,D1,f1,dplus_2,fomega_2,D2,f2);
      } else {
         fprintf(fd,"%10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e \n",
         	      a,z,omega_r,omega_m,omega_k,omega_de,hubble_parameter(a),dplus,fomega,D1,f1);
      }
  }
 
  fclose(fd);

  if (SecondOrderDplus == 1) {
     free(ScaleFactor);
     free(Growth_of_a);
     free(Growth_y2);
  }

}

