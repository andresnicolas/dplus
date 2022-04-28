/* Routines adapted from Numerical Recipies in C (Press et al. 1992) */

#include "tools.h"
#include "numrec.h"

/* Allocate an array */

double *vector(long nl, long nh)
{
   double *v;

   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) error("allocation failure in vector()");
   return v-nl+NR_END;
}

/* Free an allocated array */

void free_vector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

/* Cash-Karp-Runge-Kutta step (needed by rkqs) */

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []))
{
   int i;
   static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
                 b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
                 b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
                 b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
                 b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
                 c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
                 dc5 = -277.00/14336.0;
   double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
          dc4=c4-13525.0/55296.0,dc6=c6-0.25;
   double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
   
   ak2=vector(1,n);
   ak3=vector(1,n);
   ak4=vector(1,n);
   ak5=vector(1,n);
   ak6=vector(1,n);
   ytemp=vector(1,n);
   for (i=1;i<=n;i++)
       ytemp[i]=y[i]+b21*h*dydx[i];
   (*derivs)(x+a2*h,ytemp,ak2);
   for (i=1;i<=n;i++)
       ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
   (*derivs)(x+a3*h,ytemp,ak3);
   for (i=1;i<=n;i++)
       ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
   (*derivs)(x+a4*h,ytemp,ak4);
   for (i=1;i<=n;i++)
       ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
   (*derivs)(x+a5*h,ytemp,ak5);
   for (i=1;i<=n;i++)
       ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
   (*derivs)(x+a6*h,ytemp,ak6);
   for (i=1;i<=n;i++)
       yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
   for (i=1;i<=n;i++)
       yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
   free_vector(ytemp,1,n);
   free_vector(ak6,1,n);
   free_vector(ak5,1,n);
   free_vector(ak4,1,n);
   free_vector(ak3,1,n);
   free_vector(ak2,1,n);
}

/* Integrate one step of ODEs with accuracy monitoring (needed by odeint) */

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
   int i;
   double errmax,h,htemp,xnew,*yerr,*ytemp;
   
   yerr=vector(1,n);
   ytemp=vector(1,n);
   h=htry;
   for (;;) {
       rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
       errmax=0.0;
       for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
       errmax /= eps;
       if (errmax <= 1.0) break;
       htemp=SAFETY*h*pow(errmax,PSHRNK);
       h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
       xnew=(*x)+h;
       if (xnew == *x) error("stepsize underflow in rkqs");
   }
   if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
   else *hnext=5.0*h;
   *x += (*hdid=h);
   for (i=1;i<=n;i++) y[i]=ytemp[i];
   free_vector(ytemp,1,n);
   free_vector(yerr,1,n);
}

/* Integrate ODEs with accuracy monitoring */

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])))
{
   int nstp,i;
   double x,hnext,hdid,h;
   double *yscal,*y,*dydx;
   
   yscal=vector(1,nvar);
   y=vector(1,nvar);
   dydx=vector(1,nvar);
   x=x1;
   h=SIGN(h1,x2-x1);
   *nok = (*nbad) = 0;
   for (i=1;i<=nvar;i++) y[i]=ystart[i];
   for (nstp=1;nstp<=MAXSTP;nstp++) {
       (*derivs)(x,y,dydx);
       for (i=1;i<=nvar;i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
       if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
       (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
       if (hdid == h) ++(*nok); else ++(*nbad);
       if ((x-x2)*(x2-x1) >= 0.0) {
          for (i=1;i<=nvar;i++) ystart[i]=y[i];
   	  free_vector(dydx,1,nvar);
   	  free_vector(y,1,nvar);
   	  free_vector(yscal,1,nvar);
   	  return;
   	}
       if (fabs(hnext) <= hmin) error("Step size too small in odeint");
       h=hnext;
   }
   error("Too many steps in routine odeint");
}

/* Polynomial interpolation (needed by qromb) */

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c,*d;
   
   dif=fabs(x-xa[1]);
   c=vector(1,n);
   d=vector(1,n);
   for (i=1;i<=n;i++) {
       if ((dift=fabs(x-xa[i])) < dif) {
          ns=i;
   	  dif=dift;
       }
       c[i]=ya[i];
       d[i]=ya[i];
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
       for (i=1;i<=n-m;i++) {
           ho=xa[i]-x;
   	   hp=xa[i+m]-x;
   	   w=c[i+1]-d[i];
   	   if ((den=ho-hp) == 0.0) error("Error in routine polint");
   	   den=w/den;
   	   d[i]=hp*den;
   	   c[i]=ho*den;
   	}
   	*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   free_vector(d,1,n);
   free_vector(c,1,n);
}

/* Trapezoidal rule (needed by qromb) */

double trapzd(double (*func)(double), double a, double b, int n)
{
   double x,tnm,sum,del;
   static double s;
   int it,j;
   
   if (n == 1) {
      return (s=0.5*(b-a)*((*func)(a)+(*func)(b)));
   } else {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
          tnm=it;
   	  del=(b-a)/tnm;
   	  x=a+0.5*del;
   	  for (sum=0.0,j=1;j<=it;j++,x+=del) 
	      sum += (*func)(x);
   	  s=0.5*(s+(b-a)*sum/tnm);
   	  return s;
   }
}

/* Integrate using Romberg adaptive method */

double qromb(double (*func)(double), double a, double b, double tol)
{
   double ss,dss;
   double s[JMAXP],h[JMAXP+1];
   int j;

   h[1]=1.0;
   for (j=1;j<=JMAX;j++) {
       s[j]=trapzd(func,a,b,j);
       if (j >= KK) {
          polint(&h[j-KK],&s[j-KK],KK,0.0,&ss,&dss);
	  if (fabs(dss) <= tol*fabs(ss)) return ss;
       }
       h[j+1]=0.25*h[j];
   }
   error("Too many steps in routine qromb");
   return 0.0;
}

/* Construct a cubic spline */

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
   int i,k;
   double p,qn,sig,un,*u;
   
   u=vector(1,n-1);
   if (yp1 > 0.99e30)
      y2[1]=u[1]=0.0;
   else {
      y2[1] = -0.5;
      u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
   }
   for (i=2;i<=n-1;i++) {
       sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
       p=sig*y2[i-1]+2.0;
       y2[i]=(sig-1.0)/p;
       u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
       u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
   }
   if (ypn > 0.99e30)
      qn=un=0.0;
   else {
      qn=0.5;
      un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
   }
   y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
   for (k=n-1;k>=1;k--)
       y2[k]=y2[k]*y2[k+1]+u[k];
   free_vector(u,1,n-1);
}

/* Cubic spline interpolation */

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
  int klo,khi,k;
  double h,b,a;

  klo=1;
  khi=n;
  while (khi-klo > 1) {
     k=(khi+klo) >> 1;
     if (xa[k] > x) 
	khi=k;
     else 
	klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) error("Bad xa input to routine splint");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

