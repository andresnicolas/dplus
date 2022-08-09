
#include "io.h"
#include "cosmology.h"

int main(int argc, char **argv)
{

  double dplus,dplus_2,fit_dplus,fit_dplus_2;	
  double fomega,fomega_2,fit_fomega,fit_fomega_2;	
  double omega_m,omega_k,omega_de;

  if (argc < 2) {
      fprintf(stdout, "\n Error. Missing input file.\n");
      fprintf(stdout, "./dplus.x <input_param>\n\n");
      exit(EXIT_FAILURE);
  }

  read_inputfile(argv[1]);
  if (P.UseTab_wEoS == 1) set_dark_energy_tables();
  set_dplus_spline();

  double aini = log10(1.0/(1.0 + P.InitialRedshift));
  double aend = log10(1.0/(1.0 + P.FinalRedshift));
  double da = (aend - aini)/(double)(P.Nbins-1);

  FILE *fd = fopen(P.OutputFile,"w");

  for (int i=0; i<P.Nbins; i++) {

      double a = pow(10.0,aini + da*(double)i);
      double z = 1.0/a - 1.0;

      density_parameters(a, &omega_m, &omega_k, &omega_de);
      growth_factor(a, &dplus, &fomega);
      growth_factor_2(a, &dplus_2, &fomega_2);	
      fitting_formulas(a, &fit_dplus, &fit_dplus_2, &fit_fomega, &fit_fomega_2); 

      fprintf(fd,"%10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e %10.8e \n",
      	      a,z,omega_m,omega_k,omega_de,hubble_parameter(a),dplus,fomega,fit_dplus,fit_fomega,dplus_2,fomega_2,
              fit_dplus_2,fit_fomega_2);
  }

  if (P.UseTab_wEoS == 1) {
     free(DEtab.ScaleFactor); 
     free(DEtab.EoS);         
     free(DEtab.Factor);      
     free(DEtab.EoS_y2);      
     free(DEtab.Factor_y2);   
  }
  free(Gtab.ScaleFactor);  
  free(Gtab.Growth);       
  free(Gtab.Growth_y2);    

}

