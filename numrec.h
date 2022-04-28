#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
#define NR_END 1
#define FREE_ARG char*

#define FMAX(a,b) (a > b ? a : b)
#define FMIN(a,b) (a < b ? a : b)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []));
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])));
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd(double (*func)(double), double a, double b, int n);
double qromb(double (*func)(double), double a, double b, double tol);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
