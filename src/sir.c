#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define GAMMA       (p[parindex[0]]) // recovery rate
#define MU          (p[parindex[1]]) // death rate
#define ETA         (p[parindex[2]]) // import rate
#define BETA0       (p[parindex[3]]) // baseline transmission rate
#define BETA1       (p[parindex[4]]) // amplitude of seasonal cycle in transmission rate
#define PERIOD      (p[parindex[5]]) // period of seasonal cycle in transmission rate
#define T1          (p[parindex[6]]) // time at which pop size starts to grow
#define ALPHA       (p[parindex[7]]) // expected change in population per unit time after tbase
#define RHO         (p[parindex[8]]) // reporting probability
#define S0          (p[parindex[9]]) // initial fraction of S
#define I0          (p[parindex[10]]) // initial fraction of I
#define R0          (p[parindex[11]]) // initial fraction of R
#define N0          (p[parindex[12]]) // population size

#define SUSC      (x[stateindex[0]]) // number of susceptibles
#define INFD      (x[stateindex[1]]) // number of infectives
#define RCVD      (x[stateindex[2]]) // number of recovereds
#define POPN      (x[stateindex[3]]) // population size
#define CASE      (x[stateindex[4]]) // number of cases (accumulated per reporting period)

double _sir_rates (int j, double t, double *x, double *p,
		   int *stateindex, int *parindex, int *covindex,
		   int ncovar, double *covar) {
  double beta;
  double rate = 0.0;

  switch (j) {
  case 1: 			// birth
    rate = N0 * MU;
    if (t > T1) rate += ALPHA * (t - T1) * MU;
    break;
  case 2:			// susceptible death
    rate = MU*SUSC;
    break;
  case 3:			// infection
    beta = BETA0 * (1 + BETA1 * sinpi ( 2 * t / PERIOD));
    rate = (beta * INFD + ETA) * SUSC;
    break;
  case 4:			// infected death
    rate = MU*INFD;
    break;
  case 5:			// recovery
    rate = GAMMA*INFD;
    break;
  case 6:			// recovered death
    rate = MU*RCVD;
    break;
  default:
    error("unrecognized rate code %d",j);
    break;
  }
  return rate;
}
