
#' Get dynamic autocorrelation function.
#'
#' \code{get_dynamic_acf} estimates time-dependent autocorrelation
#' function from ensemble time series.
#'
#' Any missing values in 'x' will cause an error.
#'
#' Bandwidths affect weights in local smoothers as follows. To get the
#' local estimate corresponding to index i, the distance to each other
#' index j is calculated as (i - j) / h, where h is the
#' bandwidth. Then that distance is plugged into the kernel function
#' to obtain a weight. The weights are normalized to sum to one for
#' each index.
#'
#' The gaussian kernel is equivalent to a standard Gaussian density
#' function. The uniform kernel is an indicator function of whether
#' the distance is less than 1. Thus selecting a uniform kernel with a
#' bandwidth of 2 is equivalent to a sliding window of length 3 that
#' is centered on the focal index.
#'
#' '"local_constant"' smoothers are local means computed with the
#' kernel weights. '"local_linear"' smoothers are the fitted values of
#' local linear regressions with the kernel weights. The linear
#' smoothers avoid biases that the one-sided kernels at the ends of
#' the time series can create for the local constant smoothers.
#'
#' @param x A univariate or multivariate numeric time series object or
#' a numeric vector or matrix.
#' @param center_trend Character string giving method of calculating
#' the trend to subtract. Allowed values are '"assume_zero"',
#' '"grand_mean"', '"ensemble_means"', '"local_constant"', and
#' '"local_linear"'. Will be partially matched.
#' @param center_kernel Character string giving the kernel for any
#' local detrending. Allowed values are '"gaussian"' and '"uniform"'.
#' @param center_bandwidth Bandwith of kernel for any local detrending
#' done. A numeric value >= 1.
#' @param acf_trend Character string giving method of smoothed
#' acf. Allowed values are '"local_constant"', and
#' '"local_linear"'. Will be partially matched.
#' @param acf_kernel Character string giving the kernel for local
#' smoothing of acf. Allowed values are '"gaussian"' and '"uniform"'.
#' @param acf_bandwidth Bandwith of kernel for any local smoothing of
#' acf done. A numeric value >= 1.
#' @param acf Character string giving the type of acf to be
#' returned. Allowed values are '"correlation"' (the default), and
#' '"covariance"'. Will be partially matched.
#' @param lag Integer lag at which to calculate the acf. This lag is in terms
#' of the index of \code{x} and does not account for the frequency of
#' \code{x} if \code{x} is a time series. It should be positive.
#' @return A list with elements '"centered"' and '"acf"'. '"centered"'
#' is a list of the detrend time series, the trend subtracted, and the
#' bandwidth used in the detrending. '"acf"' is a list of the smoothed
#' estimate of the acf and the bandwidth used. If no bandwidths were
#' used they will be NULL.
#'
#' @seealso \code{\link{acf}} for regular autcorrelation estimation
#' @export
#' @examples
#'
#' # A highly autocorrelated time series
#' x <- 1:10
#' get_dynamic_acf(x, acf_bandwidth=3)
#'
#' # Plot log of acf
#' plot(log(get_dynamic_acf(x, acf_bandwidth=3)$acf$smooth))
#'
#' # Check estimates with AR1 simulations with lag-1 core 0.1
#' w <- rnorm(1000)
#' xnext <- function(xlast, w) 0.1 * xlast + w
#' x <- Reduce(xnext, x=w, init=0, accumulate=TRUE)
#' acf(x, lag.max=1, plot=FALSE)
#' head(get_dynamic_acf(x, acf_bandwidth=length(x))$acf$smooth)
#'
#' # Check detrending ability
#' x2 <- x + seq(1, 10, len=length(x))
#' ans <- get_dynamic_acf(x2, center_trend="local_linear",
#'                        center_bandwidth=length(x), acf_bandwidth=length(x))
#' head(ans$acf$smooth)
#'
#' # The simple acf estimate is inflated by the trend
#' acf(x2, lag.max=1, plot=FALSE)
#'
#'# Check ability to estimate time-dependent autocorrelation
#' xnext <- function(xlast, w) 0.8 * xlast + w
#' xhi <- Reduce(xnext, x=w, init=0, accumulate=TRUE)
#' acf(xhi, lag.max=1, plot=FALSE)
#' wt <- seq(0, 1, len=length(x))
#' xdynamic <- wt * xhi + (1 - wt) * x
#' get_dynamic_acf(xdynamic, acf_bandwidth=100)$acf$smooth
get_dynamic_acf <- function(x, center_trend="grand_mean",
                            center_kernel="gaussian",
                            center_bandwidth=NULL,
                            acf_trend="local_constant",
                            acf_kernel="uniform", acf_bandwidth=NULL,
                            acf="correlation", lag=1){
  centered <- detrend(x, trend=center_trend, kernel=center_kernel,
                      bandwidth=center_bandwidth)
  acf <- autocor(centered$x, trend=acf_trend, acf_kernel,
                 bandwidth=acf_bandwidth, cortype=acf, lag=lag)
  list(centered=centered, acf=acf)
}

detrend <- function(x, trend=c("grand_mean", "ensemble_means",
                           "local_constant", "local_linear",
                               "assume_zero"),
                    kernel=c("gaussian", "uniform"),
                    bandwidth=NULL){
  trend <- match.arg(trend)
  kernel <- match.arg(kernel)
  x <- na.fail(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric")
  rmn <- rowMeans(x)
  if (trend == "grand_mean"){
    center <- mean(rmn)
    x <- x - center
  } else if (trend == "ensemble_means"){
    center <- rmn
    x <- x - center
  } else if (trend == "local_constant"
             || trend == "local_linear"){
    samplet <- as.integer(nrow(x))
    step <- seq(1, samplet)
    data <- data.frame(step=step, rmn=rmn)
    srm <- smooth(data=data, bandwidth=bandwidth, est=trend, kernel=kernel)
    center <- srm$smooth
    x <- x - center
    bandwidth <- srm$bandwidth
  } else if (trend == "assume_zero"){
    center <- rep(0, nrow(x))
  }
  list(x=x, center=center, bandwidth=bandwidth)
}

smooth <- function(data, est, kernel="gaussian", bandwidth){
  if (!is.numeric(bandwidth) || length(bandwidth) > 1) {
    stop("argument \"bandwidth\" must be provided as a single numeric value")
  } else if (bandwidth < 1){
    stop("argument \"bandwidth\" must be >= 1")
  }
  if(kernel == "gaussian") {
    kern <- function(ind, bw=bandwidth){
      dist <- abs(data$step - ind) / bw
      w <- dnorm(dist)
      w / sum(w)
    }
  } else {
    kern <- function(ind, bw=bandwidth){
      dist <- abs(data$step - ind) / bw
      w <- dist < 1
      w / sum(w)
    }
  }
  w <- sapply(data$step, kern)
  if (est == "local_constant"){
    smooth <- colSums(w * data$rmn)
  } else {
    wlm <- function(x){
      m <- lm(rmn~step, weights=w[, x], data=data)
      fitted(m)[x]
    }
    smooth <- sapply(seq_along(data$step), wlm)
  }
  list(smooth=smooth, bandwidth=bandwidth)
}

autocor <- function(x, cortype=c("correlation", "covariance"), lag=1,
                    bandwidth=NULL, trend=c("local_constant", "local_linear"),
                    kernel=c("gaussian", "uniform")){
  trend <- match.arg(trend)
  cortype <- match.arg(cortype)
  kernel <- match.arg(kernel)
  x <- na.fail(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (lag < 0) stop("'lag' must be >= 0")

  n <- nrow(x)
  end1 <- n - lag
  start2 <- 1 + lag
  x1 <- x[1:end1, , drop=FALSE]
  x2 <- x[start2:n, , drop=FALSE]
  xx_lag <- rowMeans(x1 * x2)

  step <- seq(start2, n)
  data <- data.frame(step=step, rmn=xx_lag)
  xx_lag_sm <- smooth(data=data, bandwidth=bandwidth, kernel=kernel, est=trend)
  if (cortype == "correlation") {
    xx <- rowMeans(x1 * x1) * 0.5 + rowMeans(x2 * x2) * 0.5
    data <- data.frame(step=step, rmn= xx)
    xx_sm <- smooth(data=data, bandwidth=bandwidth, kernel=kernel, est=trend)
    ret <- list(smooth=xx_lag_sm$smooth / xx_sm$smooth,
                bandwidth=xx_sm$bandwidth)
  } else {
    ret <- xx_lag_sm
  }
  ret
}
