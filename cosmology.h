
struct eos_table {
   int    Nbins;
   double *ScaleFactor;
   double *EoS;
   double *Factor;
   double *EoS_y2;
   double *Factor_y2;
};
extern struct eos_table DEtab;

struct growth_table {
   int    Nbins;
   double *ScaleFactor;
   double *Growth;
   double *Growth_y2;
};
extern struct growth_table Gtab;

void set_dark_energy_tables(void);
void set_dplus_spline(void);
double dark_energy_eos(double a);
double dark_energy_factor_integ(double lna);
double dark_energy_factor(double a);
double hubble_parameter(double a);
void density_parameters(double a, double *omega_m, double *omega_k, double *omega_de);
void growth_factor_ode(double lna, double y[], double dydx[]);
void growth_factor(double a, double *dplus, double *fomega);
void growth_factor_2_ode(double lna, double y[], double dydx[]);
void growth_factor_2(double a, double *dplus_2, double *fomega_2);
void fitting_formulas(double a, double *D1, double *D2, double *f1, double *f2);

