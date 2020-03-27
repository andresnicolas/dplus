
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Numerical Recipies */

#define EPS 1.0e-6
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

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

extern int    kmax,kount;
extern double *xp,**yp,dxsav;

extern int    UseTab_wEOS;
extern char   OutputFile[500];
extern char   wEOSFile[500];
extern double OmegaMatter0;
extern double OmegaRadiation0;
extern double OmegaCurvature0;
extern double OmegaDarkEnergy0;
extern double Hubble0;
extern double w0EOS;
extern double waEOS;
extern double InitialRedshift;
extern double FinalRedshift;
extern int    Nbins;

extern int    NTAB;
extern double *atab,*wtab,*ftab;

