
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
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);



void read_inputfile(char *filename);
int count_lines(char *filename);
void set_dark_energy_tables(void);
void set_dplus_spline(void);
double dark_energy_eos(double a);
double dark_energy_factor_integ(double lna);
double dark_energy_factor(double a);
double hubble_parameter(double a);
void density_parameters(double a, double *omega_r, double *omega_m, double *omega_k, double *omega_de);
void growth_factor_ode(double lna, double y[], double dydx[]);
void growth_factor(double a, double *dplus, double *fomega);
void growth_factor_2_ode(double lna, double y[], double dydx[]);
void growth_factor_2(double a, double *dplus_2, double *fomega_2);
double fitting_formulas(double a, double *D1, double *D2, double *f1, double *f2);

