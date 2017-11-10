#include <R_ext/Rdynload.h>
#include "spaero_internal.h"

static const R_CallMethodDef callMethods[] = {
  {"_transition_rates_frequency_dependent", (DL_FUNC) &_transition_rates_frequency_dependent, 9},
  {"_transition_rates_density_dependent", (DL_FUNC) &_transition_rates_density_dependent, 9},
  {NULL, NULL, 0}
};

void R_init_spaero (DllInfo *info) {
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
