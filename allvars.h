
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Numerical Recipies 

#define EPS 1.0e-8
#define MAXSTP 10000
#define TINY 1.0e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define JMAX 40
#define JMAXP (JMAX+1)
#define KK 5
#define FUNC(x) ((*func)(x))
#define NR_END 1
#define FREE_ARG char*

#define FMAX(a,b) (a > b ? a : b)
#define FMIN(a,b) (a < b ? a : b)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// Parameterfile 

extern int    UseTab_wEoS;
extern int    SecondOrderDplus;
extern char   OutputFile[500];
extern char   wEoSFile[500];
extern double OmegaMatter0;
extern double OmegaRadiation0;
extern double OmegaCurvature0;
extern double OmegaDarkEnergy0;
extern double Hubble0;
extern double w0EoS;
extern double waEoS;
extern double InitialRedshift;
extern double FinalRedshift;
extern int    Nbins;

// For tabulated w(a)

extern int    SizeTable;
extern double *ScaleFactorTable;
extern double *DarkEnergyEoSTable;
extern double *DarkEnergyFactorTable;

// For order 2 growth factor

extern int    Ntab;
extern double *ScaleFactor;
extern double *Growth_of_a;
extern double *Growth_y2;
