
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
