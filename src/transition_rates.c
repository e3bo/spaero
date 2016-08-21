#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define GAMMA       (p[parindex[0]]) // recovery rate
#define MU          (p[parindex[1]]) // birth rate
#define D           (p[parindex[2]]) // death rate
#define ETA         (p[parindex[3]]) // sparking rate
#define BETA        (p[parindex[4]]) // transmission rate
#define RHO         (p[parindex[5]]) // reporting probability
#define S0          (p[parindex[6]]) // initial fraction of S
#define I0          (p[parindex[7]]) // initial fraction of I
#define R0          (p[parindex[8]]) // initial fraction of R
#define N0          (p[parindex[9]]) // initial population size

#define GAMMA_T     (covar[covindex[0]]) // recovery rate
#define MU_T        (covar[covindex[1]]) // birth rate
#define D_T         (covar[covindex[2]]) // death rate
#define ETA_T       (covar[covindex[3]]) // sparking rate
#define BETA_T      (covar[covindex[4]]) // transmission rate

#define SUSC        (x[stateindex[0]]) // number of susceptibles
#define INFD        (x[stateindex[1]]) // number of infectives
#define RCVD        (x[stateindex[2]]) // number of recovereds
#define POPN        (x[stateindex[3]]) // population size
#define CASE        (x[stateindex[4]]) // number of cases (accumulated per reporting period)

double _transition_rates_density_dependent (int j, double t, double *x, double *p,
		          int *stateindex, int *parindex, int *covindex,
		          int ncovar, double *covar) {
  double rate = 0.0;

  switch (j) {
  case 1: 			// birth
    rate = N0 * (MU + MU_T);
    break;
  case 2:			// susceptible death
    rate = SUSC * (D + D_T);
    break;
  case 3:			// infection
    rate = ((BETA + BETA_T) * INFD + (ETA + ETA_T)) * SUSC;
    break;
  case 4:			// infected death
    rate = INFD * (D + D_T);
    break;
  case 5:			// recovery
    rate = INFD * (GAMMA + GAMMA_T);
    break;
  case 6:			// recovered death
    rate = RCVD * (D + D_T);
    break;
  default:
    error("unrecognized rate code %d",j);
    break;
  }
  return rate;
}

double _transition_rates_frequency_dependent (int j, double t, double *x, double *p,
		          int *stateindex, int *parindex, int *covindex,
		          int ncovar, double *covar) {
  double rate = 0.0;

  switch (j) {
  case 1: 			// birth
    rate = N0 * (MU + MU_T);
    break;
  case 2:			// susceptible death
    rate = SUSC * (D + D_T);
    break;
  case 3:			// infection
    rate = ((BETA + BETA_T) * INFD / POPN + (ETA + ETA_T)) * SUSC;
    break;
  case 4:			// infected death
    rate = INFD * (D + D_T);
    break;
  case 5:			// recovery
    rate = INFD * (GAMMA + GAMMA_T);
    break;
  case 6:			// recovered death
    rate = RCVD * (D + D_T);
    break;
  default:
    error("unrecognized rate code %d",j);
    break;
  }
  return rate;
}
