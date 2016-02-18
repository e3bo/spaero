#' Create SIR gillespie simulator
#'
#' @useDynLib spaero
#' @export
#' @examples
#'
#' # See strong seasonal forcing
#' foo <- create_sir_simulator()
#' params <- structure(c(gamma=24, mu=0.01, eta=10e-06,
#'                       beta0=11e-06, beta1=.99,
#'                       period=1, t1=10, alpha=0,
#'                       rho=0.1, S_0=0.9974, I_0=0,
#'                       R_0=1 - .9974, N_0=1e+06))
#' out <- pomp::simulate(foo, times=seq(0, 20, by=1/26), params=params)
#' out <- as(out, "data.frame")
#'
#' opar <- par(mfrow=c(4, 1))
#' plot((S/N)~time, data=out, type="l")
#' plot(cases~time, data=out, type="l")
#' repn <- params["beta0"] * params["N_0"] / (params["mu"] + params["gamma"])
#' repn <- repn * (1 + params["beta1"] * sinpi (2 * out$time / params["period"]))
#' plot(out$time, repn, type='l')
#' plot(out$cases, repn)
#' par(opar)
#'
create_sir_simulator <- function() {
  if (!requireNamespace("pomp", quietly = TRUE)) {
    stop("The pomp package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  pomp::pomp(
     data=data.frame(time=seq(1, 10), reports=NA),
     times="time",
     t0=0,
     params=c(
       gamma=24, mu=1/70, eta=1e-4,
       beta0=9.6e-5, beta1=0, period=1,
       t1=100, alpha=100, rho=0.1,
       S_0=0.99, I_0=1e-4, R_0=0, N_0=1e5),
     rprocess=pomp::gillespie.sim(
       rate.fun="_sir_rates",
       PACKAGE="spaero",
       v=cbind(
         birth=c(1,0,0,1,0),
         sdeath=c(-1,0,0,-1,0),
         infection=c(-1,1,0,0,0),
         ideath=c(0,-1,0,-1,0),
         recovery=c(0,-1,1,0,1),
         rdeath=c(0,0,-1,-1,0)
         ),
       d=cbind(
         birth=c(1,1,1,1,0),
         sdeath=c(1,0,0,0,0),
         infection=c(1,1,1,1,0),
         ideath=c(0,1,0,0,0),
         recovery=c(0,1,0,0,0),
         rdeath=c(0,0,1,0,0)
         )
       ),
     measurement.model=reports~binom(size=cases,prob=rho),
      statenames=c("S","I","R","N","cases"),
     paramnames=c(
       "gamma","mu","eta", "beta0",
         "beta1","period","t1",
         "alpha", "rho", "S_0", "I_0", "R_0", "N_0"),
     zeronames=c("cases"),
     comp.names=c("S","I","R"),
     ic.names=c("S_0","I_0","R_0"),
     initializer=function(params, t0, comp.names, ic.names, ...) {
       x0 <- setNames(numeric(5),c("S","I","R","N","cases"))
       fracs <- params[ic.names]
       x0["N"] <- params["N_0"]
       x0[comp.names] <- round(params["N_0"] * fracs / sum(fracs))
       x0
     }
     )
}

create_pop_growth_params <- function(t_crit=300, t1=100, N_0=1e6,
                          initial_reproduction_number=0.8,
                          gamma=24, mu=1 / 70, beta1=0,
                          period=1, rho=0.1, S_0=0.99,
                          I_0=1e-6, R_0=0, eta=1/N_0){
  stopifnot(t_crit > t1)
  N_crit <- 1 / initial_reproduction_number * N_0
  alpha <- (N_crit - N_0) / (t_crit - t1)
  beta0 <- initial_reproduction_number * (gamma + mu) / N_0
  c(gamma=gamma, mu=mu, eta=eta, beta0=beta0, beta1=beta1, period=period,
    t1=t1, alpha=alpha, rho=rho, S_0=S_0, I_0=I_0, R_0=R_0, N_0=N_0)
}
