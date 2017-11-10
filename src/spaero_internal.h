#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

extern double _transition_rates_density_dependent(int j, double t, double *x, double *p, int *stateindex, int *parindex, int *covindex, int ncovar, double *covar);
extern double _transition_rates_frequency_dependent(int j, double t, double *x, double *p, int *stateindex, int *parindex, int *covindex, int ncovar, double *covar);

