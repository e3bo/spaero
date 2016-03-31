#' Create surveillance data simulator.
#'
#' \code{create_simulator} creates a pomp object that will run
#' simulations of an SIR or SIS model according to Gillespie's direct
#' method and generate simulated observations of the process.
#'
#' See the vignette "Getting started with spaero" for a description of
#' the model. The "params" argument must include all model
#' parameters. These will become the default parameters for the model
#' object. They can be overriden when the simulation is run via the
#' "params" argument of \code{pomp::simulate}. The case is the same
#' for the "times" argument.
#'
#' @return A pomp object with which simulations can be run via \code{pomp::simulate}.
#' @param times a numeric vector of increasing times at which the
#' state of the simulation will be sampled.
#' @param t0 The at which the simulation is started with state
#' variable set to the initial conditions specified via params.
#' @param process_model Character string giving the process
#' model. Allowed values are '"SIR"' and '"SIS"'.
#' @param params a named numeric vector of parameter values and
#' initial conditions.
#'
#' @seealso \code{\link[pomp]{pomp}} for documentation of pomp objects
#' @useDynLib spaero
#' @export
#' @examples
#'
#' # See strong seasonal forcing
#' foo <- create_simulator()
#' params <- c(gamma=24, mu=0.01, eta=10e-06, beta0=11e-06, beta1=.99,
#'             period=1, t1=10, alpha=0, rho=0.1, S_0=0.9974, I_0=0,
#'             R_0=1 - .9974, N_0=1e+06)
#' out <- pomp::simulate(foo, times=seq(0, 20, by=1/26), params=params)
#' out <- as(out, "data.frame")
#'
#' opar <- par(mfrow=c(4, 1))
#' plot((S/N)~time, data=out, type="l")
#' plot(cases~time, data=out, type="l")
#' repn <- params["beta0"] * out$N / (params["mu"] + params["gamma"])
#' repn <- repn * (1 + params["beta1"] *
#'                       sinpi (2 * out$time / params["period"]))
#' plot(out$time, repn, type='l')
#' plot(out$cases, repn)
#' par(opar)
#'
create_simulator <- function(times=seq(0, 9), t0=min(times),
                             process_model=c("SIR", "SIS"),
                             params=c(gamma=24, mu=1 / 70, d=1 / 70, eta=1e-5,
                                 beta=1e-4, rho=0.1, S_0=1, I_0=0, R_0=0,
                                 N_0=1e5),
                             covar=data.frame(gamma_t=c(0, 0), mu_t=c(0, 0),
                                 d_t=c(0, 0), eta_t=c(0, 0), beta_t=c(0, 0),
                                 time=c(0, 1e6))) {
  process_model <- match.arg(process_model)
  if (!requireNamespace("pomp", quietly = TRUE)) {
    stop(paste("The pomp package is needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
  data <- data.frame(time=times, reports=NA)
  d <- cbind(birth=c(1,1,1,1,0),
             sdeath=c(1,1,1,1,0),
             infection=c(1,1,1,1,0),
             ideath=c(1,1,1,1,0),
             recovery=c(1,1,1,1,0),
             rdeath=c(1,1,1,1,0))
  if (process_model == "SIR") {
    v <- cbind(birth=c(1,0,0,1,0),
               sdeath=c(-1,0,0,-1,0),
               infection=c(-1,1,0,0,0),
               ideath=c(0,-1,0,-1,0),
               recovery=c(0,-1,1,0,1),
               rdeath=c(0,0,-1,-1,0))
  } else {
    v <- cbind(birth=c(1,0,0,1,0),
               sdeath=c(-1,0,0,-1,0),
               infection=c(-1,1,0,0,0),
               ideath=c(0,-1,0,-1,0),
               recovery=c(1,-1,0,0,1),
               rdeath=c(0,0,-1,-1,0))
  }
  rprocess <- pomp::gillespie.sim(rate.fun="_transition_rates",
                                  PACKAGE="spaero", v=v, d=d)
  initializer <- function(params, t0, ...) {
    comp.names <- c("S", "I", "R")
    ic.names <- c("S_0", "I_0", "R_0")
    x0 <- setNames(numeric(5), c("S", "I", "R", "N", "cases"))
    fracs <- params[ic.names]
    x0["N"] <- params["N_0"]
    x0[comp.names] <- round(params["N_0"] * fracs / sum(fracs))
    if(params["rho"] < 0 | params["rho"] > 1) {
      stop("rho must be in [0, 1]")
    }
    pos.names <- c("gamma", "mu", "d", "eta", "beta",
                   "S_0", "I_0", "R_0", "N_0")
    if(any(params[pos.names] < 0)) {
      stop(paste("All", paste(pos.names, collapse=" "), "should be >= 0."))
    }
    x0
  }
  pomp::pomp(data=data, times="time", t0=t0, params=params, rprocess=rprocess,
             measurement.model=reports~binom(size=cases, prob=rho),
             covar=covar, statenames=c("S", "I", "R", "N", "cases"),
             paramnames=c("gamma", "mu", "d", "eta", "beta", "rho", "S_0",
                 "I_0", "R_0", "N_0"),
             covarnames=c("gamma_t", "mu_t", "d_t", "eta_t", "beta_t"),
             tcovar="time", zeronames="cases", initializer=initializer)
}

#' Create population growth parameters.
#'
#' This function calculates parameters consistent with the population
#' starting to grow at a given time point at a constant rate such
#' that the expected value of the basic reproduction number crosses
#' one from below at a later given time point. It does this by
#' setting the transmission rate, beta, and the rate of change of the
#' birth rate, alpha, accordingly. It uses a zeroth order
#' approximation that will not be accurate when the rate of change in
#' the population growth rate occurs on the same time scale as the
#' population dynamics.
#'
#' @param t_crit the time at which the the expected value of the
#' basic reproduction number should equal 1
#' @param t1 the time at which the the population size starts to grow.
#' @param N_0 the initial population size
#' @param initial_reproduction_number the initial value of the basic
#' reproduction number
#' @param gamma the recovery rate
#' @param mu the death rate
#' @param beta1 the amplitude of seasonal variation in the transmission rate
#' @param period the period of seasonal variation in the transmission rate
#' @param rho the binomial probability of a case being reported
#' @param S_0 the initial weight of susceptibles
#' @param I_0 the initial weight of infectives
#' @param R_0 the initial weight of recovereds
#' @param eta the rate of infection from outside of the population
#' @seealso \code{\link{create_simulator}}
#' @return a named vector suitable for the "params" argument to \code{create_sir_simulator}
#' @export
#' @examples
#'
#' params <- create_pop_growth_params()
#' foo <- create_simulator()
#' out <- pomp::simulate(foo, times=seq(0, 400), params=params)
#' out <- as(out, "data.frame")
#'
#' opar <- par(mfrow=c(3, 1))
#' plot(R~time, data=out, type='l')
#' plot(N~time, data=out, type='l')
#' repn <- out$N * params["beta0"] / (params["gamma"] + params["mu"])
#' plot(out$time, repn, type='l')
#' par(opar)
#'
create_pop_growth_params <- function(t_crit=300, t1=100, N_0=1e6,
                          initial_reproduction_number=0.8,
                          gamma=24, mu=1 / 70, beta1=0,
                          period=1, rho=0.1, S_0=0.99,
                          I_0=1e-6, R_0=0, eta=1 / N_0){
  stopifnot(t_crit > t1)
  stopifnot(initial_reproduction_number < 1)
  N_crit <- 1 / initial_reproduction_number * N_0
  alpha <- (N_crit - N_0) / (t_crit - t1)
  beta0 <- initial_reproduction_number * (gamma + mu) / N_0
  c(gamma=gamma, mu=mu, eta=eta, beta0=beta0, beta1=beta1, period=period,
    t1=t1, alpha=alpha, rho=rho, S_0=S_0, I_0=I_0, R_0=R_0, N_0=N_0)
}
