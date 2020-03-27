
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

void locate(double xx[], int n, double x, int *j)
{
  int ju,jm,jl;
  int ascnd;
  
  jl = 0   ;
  ju = n+1 ;
  ascnd=(xx[n-1] > xx[0]);
  while (ju-jl > 1) 
  {
    jm=(ju+jl) >> 1;
    if (x > xx[jm-1] == ascnd)
      jl = jm ;
    else
      ju = jm ;
  }
  *j = jl ;
}

// Dark energy EOS 

void set_tables(void)
{
  FILE *fd;
  int  i;
  char filename[200];

  sprintf(filename,"%s",wEOSFile);

  NTAB = count_lines(filename);
  atab = (double *) malloc(NTAB*sizeof(double));
  wtab = (double *) malloc(NTAB*sizeof(double));
  ftab = (double *) malloc(NTAB*sizeof(double));

  fd = fopen(filename,"r");

  for (i=0;i<NTAB;i++) 
      fscanf(fd,"%lf %lf \n",&atab[i],&wtab[i]);
  fclose(fd);

  for (i=0; i<NTAB; i++) 
      ftab[i] = exp(3.0*qromb(DEfactor_integ,log(atab[i]),0.0,EPS));	  
   
}

double w_a(double a)
{	
  int    indx;
  double wa;

  if (UseTab_wEOS == 0) {
     wa = w0EOS + waEOS*(1.0 - a);	  
  } else {
     locate(atab,NTAB,a,&indx);
     wa  = (wtab[indx] - wtab[indx-1])/(atab[indx] - atab[indx-1]);
     wa *= (a - atab[indx-1]);
     wa += wtab[indx-1]; 	  
  }

  return wa;
}

// Dark energy evolution factor

double DEfactor_integ(double lna)
{
  double dw,a;
 
  a = exp(lna);
  dw = 1.0 + w_a(a);

  return dw;
}

double DEfactor(double a)
{
  int    indx;	
  double f;
	
  if (UseTab_wEOS == 0) { 
     if (a == 1.0 || (w0EOS == -1.0 && waEOS == 0.0)) {
        f = 1.0;
     } else if (w0EOS != -1.0 && waEOS == 0.0) {	
        f = pow(a,-3.0*(1.0 + w0EOS));
     } else {
	f = pow(a,-3.0*(1.0 + w0EOS + waEOS + waEOS*(1.0 - a)/log(a)));     
     }
  } else {
     locate(atab,NTAB,a,&indx);
     f  = (ftab[indx] - ftab[indx-1])/(atab[indx] - atab[indx-1]);
     f *= (a - atab[indx-1]);
     f += ftab[indx-1];
  }

  return f;
}

// Hubble parameter

double hubble_a(double a)
{
  double h;

  h = OmegaMatter0/a/a/a
    + OmegaCurvature0/a/a
    + OmegaRadiation0/a/a/a/a
    + OmegaDarkEnergy0*DEfactor(a);

  h = Hubble0*sqrt(h);  

  return h;
}

// Density parameters 

void omegas_a(double a, double *omega_r, double *omega_m, double *omega_k, double *omega_de)
{
  double e;

  e = pow(Hubble0/hubble_a(a),2);
  *omega_r  = e*OmegaRadiation0/a/a/a/a;
  *omega_m  = e*OmegaMatter0/a/a/a;
  *omega_k  = e*OmegaCurvature0/a/a;
  *omega_de = e*OmegaDarkEnergy0*DEfactor(a);

}

// Linear growth factor

void growth_ode(double lna, double y[], double dydx[])
{
  double a,F1,F2;
  double omega_r,omega_m,omega_k,omega_de;

  a = exp(lna);
  omegas_a(a,&omega_r,&omega_m,&omega_k,&omega_de);

  F1 = 2.5 + 0.5*(omega_k - omega_r - 3.0*w_a(a)*omega_de);
  F2 = 2.0*omega_k + omega_r + 1.5*(1.0 - w_a(a))*omega_de;

  dydx[1] = y[2];
  dydx[2] = y[2]*F1 + y[1]*F2;
  dydx[2] = -dydx[2];

}

void growth(double a, double *dplus, double *fomega)
{
  int    N = 2;
  int    nok,nbad;
  double *ystart,hmin,hini,lna,TOL;

  lna = log(a);

  ystart = vector(1,N);
  ystart[1] = 1.0;
  ystart[2] = 0.0;

  hmin = 0.0;
  hini = EPS;

  odeint(ystart,N,-8.0,lna,EPS,hini,hmin,&nok,&nbad,growth_ode,rkqs);
  *dplus = ystart[1]*a;
  *fomega = ystart[2]/ystart[1] + 1.0;

  free_vector(ystart,1,N);
}


