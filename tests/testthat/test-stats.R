set.seed(123)

context("TDAR")

test_that("TDAR estimation works for fixed AR(1) model", {
  # Creating large sigma matrix, suppose y is length 50
  # Basis functions
  f1 <- function(n) return(rep(1, length.out = n))
  f2 <- function(n) return(seq(from = -n / 2 + 0.5, to = n / 2 - 0.5))
  
  # a values
  sigma <- 0.01
  h0 <- matrix(c(0.2, 0.02), nrow = 2)
  a <- h0[1]*f1(29) + h0[2]*f2(29)
  y <- 0
  for (i in 1:29) {
    y <- c(y[length(y)]*a[i] + rnorm(1, mean = 0, sd = sigma), y)
  }
  x <- y[-1]
  sigma_mat <- matrix(0, nrow = length(a), ncol = length(a))
  
  # Function to fill sigma_mat
  fill_matrix <- function(sigma_mat, a, sigma, t1, t2) {
    # If t1 == t2 and t1 > 0
    if (t1 == t2) {
      # Get sum_term
      i_values <- 0:(t1 - 2)
      
      # Make sure no negative values in i_values
      if (!(FALSE %in% (i_values >= 0))) {
        get_j_indices <- function(i) return(seq(from = i + 1, to = t1 - 1))
        get_aj_values <- function(j_values) sapply(j_values, function(j) a[j]^2)
        sum_term <- unlist(
          sum(
            unlist(
              lapply(lapply(lapply(i_values, get_j_indices), get_aj_values), prod)
            )
          )
        )
        
        # Get index_term
        index_term <- sigma^2 * (1 + sum_term)
      } else {
        index_term <- sigma^2
      }
      
    }
    
    # If t2 > t1 > 0
    if (t2 > t1) {
      # Get sum_term
      i_values <- 0:(t1 - 2)
      
      # Make sure there are no negatives in i_values
      if (!(FALSE %in% (i_values >= 0))) {
        get_j_indices <- function(i) return(seq(from = i + 1, to = t1 - 1))
        get_aj_values <- function(j_values) sapply(j_values, function(j) a[j]^2)
        sum_term <- unlist(
          sum(
            unlist(
              lapply(lapply(lapply(i_values, get_j_indices), get_aj_values), prod)
            )
          )
        )
        
        # Get prod_term
        prod_idx_values <- seq(from = t1, to = (t2 - 1))
        prod_term <- prod(a[prod_idx_values])
        
        # Get index_term
        index_term <- sigma^2 * (1 + sum_term) * prod_term
      } else {
        index_term <- sigma^2
      }
    }
    
    # If t1 > t2 > 0
    if (t1 > t2) {
      # Get sum_term
      i_values <- 0:(t2 - 2)
      
      # Make sure no negative values in i_values
      if (!(FALSE %in% (i_values >= 0))) {
        get_j_indices <- function(i) return(seq(from = i + 1, to = t2 - 1))
        get_aj_values <- function(j_values) sapply(j_values, function(j) a[j]^2)
        sum_term <- unlist(
          sum(
            unlist(
              lapply(lapply(lapply(i_values, get_j_indices), get_aj_values), prod)
            )
          )
        )
        
        # Get prod_term
        prod_idx_values <- seq(from = t2, to = (t1 - 1))
        prod_term <- prod(a[prod_idx_values])
        
        # Get index_term
        index_term <- sigma^2 * (1 + sum_term) * prod_term
      } else {
        index_term <- sigma^2
      }
      
    }
    
    # Return index_term
    return(index_term)
  }
  
  # Fill sigma_mat
  for (t1 in 1:length(x)) {
    for (t2 in 1:length(x)) {
      sigma_mat[t1, t2] <- fill_matrix(sigma_mat = sigma_mat, a = a, sigma = sigma, t1 = t1, t2 = t2)
    }
  }
  
  # Get M and theta
  M <- get_TDAR_M(y, f1, f2)
  V <- (M %*% sigma_mat) %*% t(M)
  theta <- get_TDAR_theta(M = M, x = x)
  z <- as.numeric((t(theta - h0) %*% solve(V)) %*% (theta - h0))
  
  # Significance test
  expect_equal((z < qchisq(p = 0.95, df = 2)), TRUE)
})

context("detrending")

test_that("Mean-based detrending works", {
  expect_equal(detrend(1:10, trend = "grand_mean")$x, matrix(1:10 - 5.5))
  expect_equal(detrend(cbind(1:10, 2:11), trend = "grand_mean")$x,
               cbind(1:10, 2:11) - 6)
  expect_equal(detrend(1:10, trend = "ensemble")$x, matrix(rep(0, 10)))
  expect_equal(detrend(cbind(1:10, 2:11), trend = "ensemble")$x,
               cbind(1:10, 2:11) - 1:10 - 0.5)
})

test_that("Kernel-based detrending works", {
  x <- runif(1:10)
  x_ind <- seq_along(x)
  expect_equal(
    detrend(x, trend = "local_constant", bandwidth = 2, kernel = "uniform")$x,
      matrix(x - stats::ksmooth(x = x_ind, y = x, "box", bandwidth = 2,
                                x.points = x_ind)$y))
  expect_equal(
    detrend(x, trend = "local_constant", bandwidth = 3, kernel = "uniform")$x,
      matrix(x - stats::ksmooth(x = x_ind, y = x, "box", bandwidth = 4,
                                x.points = x_ind)$y))
})

test_that("Skipping detrending works", {
  expect_equal(detrend(1:10, trend = "assume_zero")$x, matrix(1:10))
})

test_that("non-numeric vector input leads to errors", {
  expect_error(detrend(c(0, "0")),
               regexp = "\'x\' must be numeric")
  expect_error(get_noncentral_moments(c(10, 10), est = "local_consant",
                                      bandwidth = 1, moment_number = 0.9),
               regexp = "\'moment_number\' must be >= 1")
})

context("smoothing")

test_that("Smoothing function works as expected", {
  skip_if_not_installed("np")
  if (is.null(options("np.messages")$np.messages)) {
    options(np.messages = TRUE)
  }
  if (is.null(options("np.tree")$np.tree)) {
    options(np.tree = FALSE)
  }
  np_smooth <- function(data, est, bandwidth, kernel = "gaussian"){
    is.constant <- grepl("constant", est)
    rt <- ifelse(is.constant, "lc", "ll")
    bw <- np::npregbw(formula = rmn ~ step, bws = bandwidth,
                      regtype = rt, ckertype = kernel,
                      bandwidth.compute = FALSE, data = data,
                      na.action = na.fail)
    mod <- np::npreg(bw)
    list(smooth = fitted(mod), bandwidth = bw$bw)
  }
  data <- data.frame(step = 1:10, rmn = 1:10)
  expect_equal(smooth(data, est = "local_constant", bandwidth = 2)$smooth,
               np_smooth(data, est = "local_constant",  bandwidth = 2)$smooth)
  expect_equal(smooth(data, est = "local_linear", bandwidth = 2)$smooth,
               np_smooth(data, est = "local_linear",  bandwidth = 2)$smooth,
               check.names = FALSE)

  expect_equal(smooth(data, est = "local_constant", kernel = "uniform",
                      bandwidth = 3)$smooth,
               np_smooth(data, est = "local_constant", kernel = "uniform",
                         bandwidth = 3)$smooth)
  expect_equal(smooth(data, est = "local_linear", bandwidth = 3)$smooth,
               np_smooth(data, est = "local_linear", bandwidth = 3)$smooth,
               check.names = FALSE)

  data2 <- data.frame(step = 20:100, rmn = rnorm(81))
  expect_equal(smooth(data2, est = "local_constant", bandwidth = 5)$smooth,
               np_smooth(data2, est = "local_constant", bandwidth = 5)$smooth)
  expect_equal(smooth(data2, est = "local_linear", bandwidth = 5,
                      kernel = "uniform")$smooth,
               np_smooth(data2, est = "local_linear", bandwidth = 5,
                         kernel = "uniform")$smooth,
               check.names = FALSE)
})

context("moment function")

test_that("argument checking works", {
            expect_error(get_noncentral_moments(c(0, "0")),
                         regexp = "\'x\' must be numeric")

          })


context("smoothing arguments")

test_that("invalid bandwidths lead to errors", {
  data <- data.frame(step = 1:10, rmn = 1:10)
  expect_error(smooth(data),
    regexp = "argument \"bandwidth\" is missing, with no default")
  expect_error(smooth(data, bandwidth = c(1:20)),
               regexp = paste("argument \"bandwidth\" must be provided as a",
                   "single numeric value"))
  expect_error(smooth(data, bandwidth = "hmm"),
               regexp = paste("argument \"bandwidth\" must be provided as a",
                   "single numeric value"))
  expect_error(smooth(data, bandwidth = 0.5),
    regexp = "argument \"bandwidth\" must be >= 1")
})

context("autocor")

test_that("invalid arguments lead to errors", {
  expect_error(autocor(NA))
  expect_error(autocor(letters), regexp = "'x' must be numeric")
  expect_error(autocor(1:10, lag = -1), regexp = "'lag' must be >= 0")
})

test_that("large bandwidth autocor estimates agree with acf", {
  skip_on_cran()

  w <- rnorm(1000)
  xnext <- function(xlast, w) 0.1 * xlast + w
  x <- Reduce(xnext, x = w, init = 0, accumulate = TRUE)

  stats_est <- stats::acf(x, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
  spaero_est <- autocor(x, est = "local_constant", kernel = "gaussian",
                        bandwidth = length(x) * 10)$smooth[2]
  expect_equal(stats_est, spaero_est, tolerance = 0.05)

  stats_est <- stats::acf(x, lag.max = 1, plot = FALSE,
                          type = "covariance")$acf[2, 1, 1]
  spaero_est <- autocor(x, bandwidth = length(x) * 10, est = "local_constant",
                        kernel = "gaussian", cortype = "covariance")$smooth[2]
  expect_equal(stats_est, spaero_est, tolerance = 0.05)

  xnext <- function(xlast, w) 0.9 * xlast + w
  x <- Reduce(xnext, x = w, init = 0, accumulate = TRUE)

  stats_est <- stats::acf(x, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
  spaero_est <- autocor(x, bandwidth = length(x) * 10, est = "local_constant",
                        kernel = "gaussian")$smooth[2]
  expect_equal(stats_est, spaero_est, tolerance = 0.05)
})

context("expected use of get_stats")

test_that("lag-0 results are sensible", {
  bw <- 10
  n <- 100
  x <- rnorm(n)
  sp <- get_stats(x, center_kernel = "uniform", center_trend = "local_constant",
                  center_bandwidth = bw, stat_bandwidth = bw,
                  stat_kernel = "uniform", lag = 0)
  expect_equal(sp$stats$autocovariance, sp$stats$variance)
  expect_equal(sp$stats$autocorrelation, rep(1, n))
})

test_that(paste("estimate of ensemble stats consistent",
                "in case of time-dependent AR(1) model"), {
  make_time_dependent_updater <- function(f) {
    t <- 0
    updater <- function(xlast, w) {
      phi <- f(t)
      t <<- t + 1
      phi * xlast + w
    }
    return(updater)
  }
  phi_t <- function(t) {
    lambda <- min(-1  + 0.01 * t, 0)
    exp (lambda)
  }
  nensemble <- 3e3
  nobs <- 90
  x <- matrix(nrow = nobs, ncol = nensemble)
  for (i in seq(1, nensemble)){
    updater <- make_time_dependent_updater(f = phi_t)
    w <- rnorm(nobs - 1)
    init <- rnorm(n = 1, sd = sqrt(1.16))
    x[, i] <- Reduce(updater, x = w, init = init, accumulate = TRUE)
  }

  est <- get_stats(x, center_trend = "assume_zero", stat_bandwidth = 3,
                   stat_trend = "local_linear")
  lambda_ests <- log(est$stats$autocorrelation[-1])
  lambda_known <- seq(from = -1, by = 0.01, len = nobs)
  var_known <-  1 / (1 - exp(2 * lambda_known))
  error_in_mean <- est$stats$mean

  expect_equal(est$centered$x + est$centered$center, x)
  expect_equal(mean(error_in_mean), 0)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.01)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm = TRUE), 0.001)
  expect_lt(mean(acov_error ^ 2, na.rm = TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt(mean(decay_time_error ^ 2), 0.5)

  expect_lt(mean( (est$stats$variance - var_known) ^ 2), 0.1)
  expect_lt(mean(est$stats$skewness ^ 2), 0.01)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.01)

  trend <- sin(2 * pi * (1:nobs) / nobs) + 2
  xx <- x + trend
  est <- get_stats(xx, center_trend = "local_constant", center_bandwidth = 3,
                   stat_bandwidth = 3)
  lambda_ests <- log(est$stats$autocorrelation[-1])

  error_in_mean <- est$stats$mean - trend
  expect_equal(est$centered$x + est$centered$center, xx)
  expect_lt(mean(error_in_mean ^ 2), 0.01)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.01)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm = TRUE), 0.001)
  expect_lt(mean(acov_error ^ 2, na.rm = TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt(mean(decay_time_error ^ 2), 0.5)

  expect_lt(mean( (est$stats$variance - var_known) ^ 2), 0.1)
  error_in_cv <- est$stats$coefficient_of_variation -  sqrt(var_known) / trend
  expect_lt(sqrt(mean( (error_in_cv) ^ 2)), 0.05)
  error_in_id <- est$stats$index_of_dispersion -  var_known / trend
  expect_lt(sqrt(mean( (error_in_id) ^ 2)), 0.2)
  expect_lt(mean(est$stats$skewness ^ 2), 0.01)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.01)
})

test_that(paste("estimate of stats consistent",
                "in case of stationary AR(1) model"), {
  skip_on_cran()

  nobs <- 4e3
  w <- rnorm(nobs - 1)
  phi <- 0.1
  xnext <- function(xlast, w) phi * xlast + w
  x <- Reduce(xnext, x = w, init = 0, accumulate = TRUE)

  est <- get_stats(x, center_trend = "assume_zero", stat_bandwidth = nobs,
                   stat_kernel = "uniform", stat_trend = "local_linear")
  lambda_ests <- log(est$stats$autocorrelation[-1])
  lambda_known <- rep(log(phi), nobs)
  var_known <-  1 / (1 - exp(2 * lambda_known))
  error_in_mean <- est$stats$mean

  expect_equal(est$centered$x[, 1] + est$centered$center, x)
  expect_equal(mean(error_in_mean), 0)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.1)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm = TRUE), 0.002)
  expect_lt(mean(acov_error ^ 2, na.rm = TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt(mean(decay_time_error ^ 2), 0.5)

  expect_lt(mean( (est$stats$variance - var_known) ^ 2), 0.1)
  expect_lt(mean(est$stats$skewness ^ 2), 0.02)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.05)

  trend <- sin(2 * pi * (1:nobs) / nobs) + 2
  xx <- x + trend
  est <- get_stats(xx, center_trend = "local_constant",
                   center_bandwidth = nobs / 16,
                   stat_bandwidth = nobs)
  lambda_ests <- log(est$stats$autocorrelation[-1])

  error_in_mean <- est$stats$mean - trend
  expect_equal(est$centered$x[, 1] + est$centered$center, xx)
  expect_lt(mean(error_in_mean ^ 2), 0.05)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.2)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm = TRUE), 0.002)
  expect_lt(mean(acov_error ^ 2, na.rm = TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt(mean(decay_time_error ^ 2), 0.5)

  expect_lt(mean( (est$stats$variance - var_known) ^ 2), 0.1)
  error_in_cv <- est$stats$coefficient_of_variation -  sqrt(var_known) / trend
  expect_lt(sqrt(mean( (error_in_cv) ^ 2)), 0.2)
  error_in_id <- est$stats$index_of_dispersion -  var_known / trend
  expect_lt(sqrt(mean( (error_in_id) ^ 2)), 0.2)
  expect_lt(mean(est$stats$skewness ^ 2), 0.01)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.02)
})

test_that(paste("Estimate of stats consistent with other methods",
                "in case of moving window estimates in",
                "nonstationary AR(1) model"), {
  skip_if_not_installed("pomp")
  skip_if_not_installed("earlywarnings")
  skip_on_cran()

  params <- c(gamma = 24, mu = 0.014, d = 0.014, eta = 1e-4, beta_par = 0,
              rho = 0.9, S_0 = 1, I_0 = 0, R_0 = 0, N_0 = 1e5, p = 0)
  covar <- data.frame(gamma_t = c(0, 0), mu_t = c(0, 0), d_t = c(0, 0),
                      eta_t = c(0, 0), beta_par_t = c(0, 24e-5), p_t = c(0, 0),
                      time = c(0, 300))
  times <- seq(0, 200, by = 1 / 12)

  sim <- create_simulator(params = params, times = times, covar = covar)
  so <- pomp::simulate(sim, as.data.frame = TRUE, seed = 272)

  bw <- 720
  n <- nrow(so)
  sp <- get_stats(diff(so[, "reports"]), center_kernel = "uniform",
                  center_trend = "local_constant", center_bandwidth = bw,
                  stat_bandwidth = bw, stat_kernel = "uniform")
  mw <- 2 * (bw - 1) + 1
  spbw <- get_stats(diff(so[, "reports"]), center_kernel = "uniform",
                    center_trend = "local_constant", center_bandwidth = mw,
                    stat_bandwidth = mw, stat_kernel = "uniform",
                    backward_only = TRUE)
  ew <- earlywarnings::generic_ews(diff(so[, "reports"]),
                                   winsize = mw / n * 100,
                                   detrending = "no")
  spm <- lapply(sp$stats, function(x) x[(bw):(n - bw)])
  spbwm <- lapply(spbw$stats, function(x) x[-seq_len(mw - 1)])
  expect_equal(spm$mean, spbwm$mean)
  expect_equal(as.data.frame(spm), as.data.frame(spbwm), tolerance = 0.05)

  ## windows sizes and hence output lengths mismatch due to rounding
  ## inside generic_ews
  ew <- ew[ew$timeindex > mw - 1, ]

  expect_equal(ew$acf1, spm$autocorrelation, tolerance = 0.01)
  expect_equal(ew$sd, sqrt(spm$variance), tolerance = 0.01)
  expect_equal(ew$kurt, spm$kurtosis, tolerance = 0.01)
  expect_equal(ew$sk, abs(spm$skewness), tolerance = 0.06)
})
