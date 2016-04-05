
context("Gillespie direct method simulator")

test_that(paste("Mean and stddev of stationary model over time",
                "consistent with ensemble mean and stdev of dizzy",
                "progam's implementation"), {

  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=24e-2,
              rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e2)
  covar <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0), eta_t=c(0, 0),
                      beta_t=c(0, 0), time=c(0, 1e6))
  times <- seq(0, 1e6, by=1)

  sim <- create_simulator(params=params, times=times, covar=covar,
                          process_model="SIS")
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=200)
  expect_lt(abs(mean(so[, "I"]) - 0.014375), 0.5)
  expect_lt(abs(mean(so[, "S"]) - 99.9888), 1)
  expect_lt(abs(sd(so[, "I"]) - 0.5130), 0.25)
  expect_lt(abs(sd(so[, "S"]) - 9.9963), 0.1)
})


test_that(paste("Means and final stddev of time-dependent model",
                "consistent with ensemble mean and stdev of dizzy",
                "progam's implementation"), {

  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=0e-2,
              rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e2)
  covar <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0), eta_t=c(0, 0),
                      beta_t=c(0, 3 * 24e-2), time=c(0, 60))
  times <- seq(0, 50, len=100)

  sim <- create_simulator(params=params, times=times, covar=covar,
                          process_model="SIS")
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=200, nsim=1000)
  ens_infected <- unstack(so, I~sim)
  ens_susceptible <- unstack(so, S~sim)
  dzout <- read.csv(file.path("dizzy", "out-linear-trend.csv"), nrows=100)
  expect_lt(sqrt(mean((dzout$I - rowMeans(ens_infected)) ^ 2)), 1)
  expect_lt(sqrt(mean((dzout$S - rowMeans(ens_susceptible)) ^ 2)), 1)
  dzout_fluc <- read.csv(file.path("dizzy", "out-linear-trend.csv"), skip=101,
                         header=FALSE)
  tfsds <- dzout_fluc[, 2]
  names(tfsds) <- dzout_fluc[, 1]
  expect_equal(sd(ens_infected[100, ]), tfsds["I"],
               check.attributes=FALSE, tol=0.1)
  expect_equal(sd(ens_susceptible[100, ]), tfsds["S"],
               check.attributes=FALSE, tol=0.1)
})
