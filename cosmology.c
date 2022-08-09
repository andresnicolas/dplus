
#include "tools.h"
#include "io.h"
#include "cosmology.h"
#include "numrec.h"

struct eos_table DEtab;
struct growth_table Gtab;


/* Dark energy equation of state */ 

void set_dark_energy_tables(void)
{

  DEtab.Nbins = count_lines(P.wEoSFile);
  DEtab.ScaleFactor = (double *) malloc(DEtab.Nbins*sizeof(double));
  DEtab.EoS = (double *) malloc(DEtab.Nbins*sizeof(double));
  DEtab.Factor = (double *) malloc(DEtab.Nbins*sizeof(double));
  DEtab.EoS_y2 = (double *) malloc(DEtab.Nbins*sizeof(double));
  DEtab.Factor_y2 = (double *) malloc(DEtab.Nbins*sizeof(double));

  FILE *fd = fopen(P.wEoSFile,"r");

  for (int i=0;i<DEtab.Nbins;i++) 
      fscanf(fd,"%lf %lf ",&DEtab.ScaleFactor[i],&DEtab.EoS[i]);
  fclose(fd);

  for (int i=0; i<DEtab.Nbins; i++) 
      DEtab.Factor[i] = exp(3.0*qromb(dark_energy_factor_integ,log(DEtab.ScaleFactor[i]),0.0,EPS));	  
 
  spline(DEtab.ScaleFactor, DEtab.EoS, DEtab.Nbins, 1.0e30, 1.0e30, DEtab.EoS_y2);   
  spline(DEtab.ScaleFactor, DEtab.Factor, DEtab.Nbins, 1.0e30, 1.0e30, DEtab.Factor_y2);   

}

double dark_energy_eos(double a)
{	
  double eos;

  if (P.UseTab_wEoS == 0) 
     eos = P.w0EoS + P.waEoS*(1.0 - a);	  
  else 
     splint(DEtab.ScaleFactor, DEtab.EoS, DEtab.EoS_y2, DEtab.Nbins, a, &eos);
  
  return eos;
}

/* Dark energy evolution factor */

double dark_energy_factor_integ(double lna)
{
  double a = exp(lna);
  return 1.0 + dark_energy_eos(a);
}

double dark_energy_factor(double a)
{
  double f;
	
  if (P.UseTab_wEoS == 0) {
     if (a == 1.0 || (P.w0EoS == -1.0 && P.waEoS == 0.0)) {
        f = 1.0;
     } else if (P.w0EoS != -1.0 && P.waEoS == 0.0) {	
        f = pow(a,-3.0*(1.0 + P.w0EoS));
     } else {
	f = pow(a,-3.0*(1.0 + P.w0EoS + P.waEoS + P.waEoS*(1.0 - a)/log(a)));     
     }
  } else 
     splint(DEtab.ScaleFactor, DEtab.Factor, DEtab.Factor_y2, DEtab.Nbins, a, &f);
  
  return f;
}

/* Hubble parameter */

double hubble_parameter(double a)
{

  double h = P.OmegaMatter / pow(a,3)
           + (1.0 - P.OmegaMatter - P.OmegaDarkEnergy) / pow(a,2)
           + P.OmegaDarkEnergy * dark_energy_factor(a);

  h = P.Hubble*sqrt(h);  

  return h;
}

/* Density parameters */ 

void density_parameters(double a, double *omega_m, double *omega_k, double *omega_de)
{
  double evol = pow(P.Hubble/hubble_parameter(a),2);

  *omega_m  = evol * P.OmegaMatter / pow(a,3);
  *omega_k  = evol * (1.0 - P.OmegaMatter - P.OmegaDarkEnergy) / pow(a,2);
  *omega_de = evol * P.OmegaDarkEnergy * dark_energy_factor(a);
}

/* Growth factor D(a) and f = dlogD/dloga - 1st order */

void growth_factor_ode(double lna, double y[], double dydx[])
{
  double a = exp(lna);
  double omega_m,omega_k,omega_de;

  density_parameters(a,&omega_m,&omega_k,&omega_de);
 
  double F1 = 2.5 + 0.5*(omega_k - 3.0*dark_energy_eos(a)*omega_de);
  double F2 = 2.0*omega_k + 1.5*(1.0 - dark_energy_eos(a))*omega_de;

  dydx[1] = y[2];
  dydx[2] = - y[2]*F1 - y[1]*F2;

}

void growth_factor(double a, double *dplus, double *fomega)
{
  int    N = 2;
  int    nok,nbad;
  double *ystart;

  ystart = vector(1,N);
  ystart[1] = 1.0;
  ystart[2] = 0.0;

  double aini = 1.0e-3;
  double hmin = 0.0;
  double hini = EPS;

  odeint(ystart,N,log(aini),log(a),EPS,hini,hmin,&nok,&nbad,growth_factor_ode,rkqs);
  *dplus = ystart[1]*a;
  *fomega = ystart[2]/ystart[1] + 1.0;

  free_vector(ystart,1,N);
}

/* Growth factor D2(a) and f2 = dlogD2/dloga - 2nd order */

void set_dplus_spline(void) 
{
  Gtab.Nbins = P.Nbins*10;
  double aini = log10(1.0e-3);
  double aend = log10(1.0/(1.0 + P.FinalRedshift));
  double da = (aend - aini)/(double)Gtab.Nbins;
  double fomega;
  
  Gtab.ScaleFactor = (double *) malloc(Gtab.Nbins*sizeof(double));
  Gtab.Growth = (double *) malloc(Gtab.Nbins*sizeof(double));
  Gtab.Growth_y2 = (double *) malloc(Gtab.Nbins*sizeof(double));

  for (int i=0; i<Gtab.Nbins; i++) {
      Gtab.ScaleFactor[i] = pow(10.0,aini + da*(double)(i+1));
      growth_factor(Gtab.ScaleFactor[i], &Gtab.Growth[i], &fomega);	  
  }

  spline(Gtab.ScaleFactor, Gtab.Growth, Gtab.Nbins, 1.0e30, 1.0e30, Gtab.Growth_y2);   

}

void growth_factor_2_ode(double lna, double y[], double dydx[])
{
  double a = exp(lna);
  double omega_m,omega_k,omega_de,dplus;

  density_parameters(a,&omega_m,&omega_k,&omega_de);
  
  splint(Gtab.ScaleFactor, Gtab.Growth, Gtab.Growth_y2, Gtab.Nbins, a, &dplus);

  double F1 = 2.5 + 0.5*(omega_k - 3.0*dark_energy_eos(a)*omega_de);
  double F2 = 2.0*omega_k + 1.5*(1.0 - dark_energy_eos(a))*omega_de;
  double F3 = 1.5*omega_m*pow(dplus,2)/a;

  dydx[1] = y[2];
  dydx[2] = - y[2]*F1 - y[1]*F2 - F3;

}

void growth_factor_2(double a, double *dplus_2, double *fomega_2)
{
  int    N = 2;
  int    nok,nbad;
  double *ystart;
  double aini = 1.0e-3;
  double hmin = 0.0;
  double hini = EPS;

  ystart = vector(1,N);
  ystart[1] = -(3.0/7.0)*aini;
  ystart[2] = (70./9.0)*aini;

  odeint(ystart,N,log(aini),log(a),EPS,hini,hmin,&nok,&nbad,growth_factor_2_ode,rkqs);
  *dplus_2 = ystart[1]*a;
  *fomega_2 = ystart[2]/ystart[1] + 1.0;

  free_vector(ystart,1,N);
}

/* Fitting formulas (Bernardeau et al. 2002, Physics Reports, 367, 1).
   These formulas are truly valid for LCDM cosmologies. */

void fitting_formulas(double a, double *D1, double *D2, double *f1, double *f2)
{
 
  double omega_m, omega_k, omega_de;

  density_parameters(a, &omega_m, &omega_k, &omega_de);

  // Growth factor

  // First order  
  *D1 = 2.5 * a * omega_m / (pow(omega_m,4.0/7.0) - omega_de + (1.0 + omega_m/2.0)*(1.0 + omega_de/70.0));
  // Second order 
  *D2 = -(3.0/7.0) * pow((*D1),2) * pow(omega_m,-1.0/143.0);

  // f = dlnD/dlna

  // First order 
  *f1 = pow(omega_m,5.0/9.0);
  // Second order 
  *f2 = 2.0*pow(omega_m,6.0/11.0);

}

