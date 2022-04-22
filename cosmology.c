
#include "allvars.h"
#include "proto.h"

int count_lines(char *filename) 
{
    FILE *fp; 
    int  count = 0;
    char string[256];

    fp = fopen(filename, "r");
    if (fp == NULL) {
       fprintf(stdout,"Error. %s doesn't exist \n",filename);
       exit(EXIT_FAILURE);
    }
    
    while (fgets(string,256,fp)) count++;
    fclose(fp);

    return count;
}

// Dark energy EoS 

void set_dark_energy_tables(void)
{
  FILE *fd;
  int  i;
  char filename[200];

  sprintf(filename,"%s",wEoSFile);

  SizeTable = count_lines(filename);
  ScaleFactorTable = (double *) malloc(SizeTable*sizeof(double));
  DarkEnergyEoSTable = (double *) malloc(SizeTable*sizeof(double));
  DarkEnergyFactorTable = (double *) malloc(SizeTable*sizeof(double));

  fd = fopen(filename,"r");

  for (i=0;i<SizeTable;i++) 
      fscanf(fd,"%lf %lf ",&ScaleFactorTable[i],&DarkEnergyEoSTable[i]);
  fclose(fd);

  for (i=0; i<SizeTable; i++) 
      DarkEnergyFactorTable[i] = exp(3.0*qromb(dark_energy_factor_integ,log(ScaleFactorTable[i]),0.0,EPS));	  
   
}

double dark_energy_eos(double a)
{	
  int    indx;
  double y2,y1,x2,x1,eos;

  if (UseTab_wEoS == 0) {

     eos = w0EoS + waEoS*(1.0 - a);	  

  } else {

     locate(ScaleFactorTable, SizeTable, a, &indx);

     y2 = DarkEnergyEoSTable[indx];
     y1 = DarkEnergyEoSTable[indx-1];
     x2 = ScaleFactorTable[indx];
     x1 = ScaleFactorTable[indx-1];

     eos = (y2 - y1) / (x2 - x1) * (a - x1) + y1;
  
  }

  return eos;
}

// Dark energy evolution factor

double dark_energy_factor_integ(double lna)
{
  double a = exp(lna);
  return 1.0 + dark_energy_eos(a);
}

double dark_energy_factor(double a)
{
  int    indx;	
  double f,y2,y1,x2,x1;
	
  if (UseTab_wEoS == 0) {

     if (a == 1.0 || (w0EoS == -1.0 && waEoS == 0.0)) {
        f = 1.0;
     } else if (w0EoS != -1.0 && waEoS == 0.0) {	
        f = pow(a,-3.0*(1.0 + w0EoS));
     } else {
	f = pow(a,-3.0*(1.0 + w0EoS + waEoS + waEoS*(1.0 - a)/log(a)));     
     }

  } else {

     locate(ScaleFactorTable, SizeTable, a, &indx);

     y2 = DarkEnergyFactorTable[indx];
     y1 = DarkEnergyFactorTable[indx-1];
     x2 = ScaleFactorTable[indx];
     x1 = ScaleFactorTable[indx-1];

     f = (y2 - y1) / (x2 - x1) * (a - x1) + y1;

  }

  return f;
}

// Hubble parameter

double hubble_parameter(double a)
{
  double h;

  h = OmegaMatter0 / pow(a,3)
    + OmegaCurvature0 / pow(a,2)
    + OmegaRadiation0 / pow(a,4)
    + OmegaDarkEnergy0 * dark_energy_factor(a);

  h = Hubble0*sqrt(h);  

  return h;
}

// Density parameters 

void density_parameters(double a, double *omega_r, double *omega_m, double *omega_k, double *omega_de)
{
  double evol = pow(Hubble0/hubble_parameter(a),2);

  *omega_r  = evol * OmegaRadiation0 / pow(a,4);;
  *omega_m  = evol * OmegaMatter0 / pow(a,3);
  *omega_k  = evol * OmegaCurvature0 / pow(a,2);
  *omega_de = evol * OmegaDarkEnergy0 * dark_energy_factor(a);
}

// Growth factor - 1st order

void growth_factor_ode(double lna, double y[], double dydx[])
{
  double a,F1,F2;
  double omega_r,omega_m,omega_k,omega_de;

  a = exp(lna);
  density_parameters(a,&omega_r,&omega_m,&omega_k,&omega_de);

  F1 = 2.5 + 0.5*(omega_k - omega_r - 3.0*dark_energy_eos(a)*omega_de);
  F2 = 2.0*omega_k + omega_r + 1.5*(1.0 - dark_energy_eos(a))*omega_de;

  dydx[1] = y[2];
  dydx[2] = - y[2]*F1 - y[1]*F2;

}

void growth_factor(double a, double *dplus, double *fomega)
{
  int    N = 2;
  int    nok,nbad;
  double *ystart,hmin,hini,aini;

  ystart = vector(1,N);
  ystart[1] = 1.0;
  ystart[2] = 0.0;

  aini = 1.0e-5;
  hmin = 0.0;
  hini = EPS;

  odeint(ystart,N,log(aini),log(a),EPS,hini,hmin,&nok,&nbad,growth_factor_ode,rkqs);
  *dplus = ystart[1]*a;
  *fomega = ystart[2]/ystart[1] + 1.0;

  free_vector(ystart,1,N);
}

// Growth factor - 2nd order

void set_dplus_spline(void) 
{
  double aini,aend,da;
  int    i;
  double fomega;

  Ntab = Nbins*10;
  aini = -6.0;
  aend = log10(1.0/(1.0 + FinalRedshift));
  da = (aend - aini)/(double)Ntab;
  
  ScaleFactor = (double *) malloc(Ntab*sizeof(double));
  Growth_of_a = (double *) malloc(Ntab*sizeof(double));
  Growth_y2 = (double *) malloc(Ntab*sizeof(double));

  for (i=0; i<Ntab; i++) {
      ScaleFactor[i] = pow(10.0,aini + da*(double)(i+1));
      growth_factor(ScaleFactor[i], &Growth_of_a[i], &fomega);	  
  }

  spline(ScaleFactor, Growth_of_a, Ntab, 1.0e30, 1.0e30, Growth_y2);   

}

void growth_factor_2_ode(double lna, double y[], double dydx[])
{
  int    indx;	
  double dplus,a,F1,F2,F3,fomega;
  double omega_r,omega_m,omega_k,omega_de;

  a = exp(lna);
  density_parameters(a,&omega_r,&omega_m,&omega_k,&omega_de);
  
  splint(ScaleFactor, Growth_of_a, Growth_y2, Ntab, a, &dplus);

  F1 = 2.5 + 0.5*(omega_k - omega_r - 3.0*dark_energy_eos(a)*omega_de);
  F2 = 2.0*omega_k + omega_r + 1.5*(1.0 - dark_energy_eos(a))*omega_de;
  F3 = 1.5*omega_m*pow(dplus,2)/a;

  dydx[1] = y[2];
  dydx[2] = - y[2]*F1 - y[1]*F2 - F3;

}

void growth_factor_2(double a, double *dplus_2, double *fomega_2)
{
  int    N = 2;
  int    nok,nbad;
  double *ystart,hmin,hini,aini;

  aini = 1.0e-5;

  ystart = vector(1,N);
  ystart[1] = -(3.0/7.0)*aini;
  ystart[2] = (70./9.0)*aini;

  hmin = 0.0;
  hini = EPS;

  odeint(ystart,N,log(aini),log(a),EPS,hini,hmin,&nok,&nbad,growth_factor_2_ode,rkqs);
  *dplus_2 = ystart[1]*a;
  *fomega_2 = ystart[2]/ystart[1] + 1.0;

  free_vector(ystart,1,N);
}

// Fitting formulas (Bernardeau et al. 2002, Physics Reports, 367, 1).
// WARNING: These formulas are truly valid for LCDM cosmologies.

double fitting_formulas(double a, double *D1, double *D2, double *f1, double *f2)
{
 
  double omega_r, omega_m, omega_k, omega_de;

  density_parameters(a, &omega_r, &omega_m, &omega_k, &omega_de);

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

