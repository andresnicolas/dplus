
/* Numerical Recipies routines */
void nrerror(char error_text[]);
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


void read_inputfile(char *filename);
void locate(double xx[], int n, double x, int *j);
int count_lines(char *filename);
void set_tables(void);
double w_a(double a);
double DEfactor_integ(double lna);
double DEfactor(double a);
double hubble_a(double a);
void omegas_a(double a, double *omega_r, double *omega_m, double *omega_k, double *omega_de);
void growth_ode(double lna, double y[], double dydx[]);
void growth(double a, double *dplus, double *fomega);
