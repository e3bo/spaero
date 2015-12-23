
#' Estimate time-dependent autocorrelation function from ensemble time series
#'
#' @param x A univariate or multivariate numeric time series object or
#' a numeric vector or matrix.
#' @param center_trend Character string giving method of calculating
#' the trend to subtract. Allowed values are '"grand_mean"',
#' '"ensemble_means"', '"local_constant"', and '"local_linear"'. Will
#' be partially matched.
#' @param center_kernel Character string giving the kernel for any
#' local detrending. Allowed values are '"gaussian"' and '"uniform"'.
#' @param center_bandwidth Bandwith of kernel for any local detrending
#' done. If not supplied, a bandwidth will be selected by
#' cross-validation of least squared error.
#' @param acf_trend Character string giving method of smoothed
#' acf. Allowed values are '"local_constant"', and
#' '"local_linear"'. Will be partially matched.
#' @param acf_kernel Character string giving the kernel for local
#' smoothing of acf. Allowed values are '"gaussian"' and '"uniform"'.
#' @param acf_bandwidth Bandwith of kernel for any local smoothing of
#' acf done. If not supplied, a bandwidth will be selected by
#' cross-validation of least squared error.
#' @param acf Character string giving the type of acf to be
#' returned. Allowed values are '"correlation"' (the default), and
#' '"covariance"'. Will be partially matched.
#' @param lag Integer lag at which to calculate the acf. This lag is in terms
#' of the index of \code{x} and does not account for the frequency of
#' \code{x} if \code{x} is a time series. It should be positive.
#' @param nmulti Integer giving the number of starting pionts in
#' search for bandwidths if any are selected by cross validation.
#'
#' @export
#'
#' @details Any missing values in 'x' will cause an error. bandwidths
get_dynamic_acf <- function(x, center_trend="grand_mean",
                            center_kernel="gaussian",
                            center_bandwidth=NULL,
                            acf_trend="local_constant",
                            acf_kernel="uniform", acf_bandwidth=NULL,
                            acf="correlation", lag=1, nmulti=5){
  centered <- detrend(x, trend=center_trend, kernel=center_kernel,
                      bandwidth=center_bandwidth, nmulti=nmulti)
  acf <- autocor(centered$x, trend=acf_trend, acf_kernel,
                 bandwidth=acf_bandwidth, cortype=acf, lag=lag)
  list(centered=centered, acf=acf)
}

detrend <- function(x, trend=c("grand_mean", "ensemble_means",
                           "local_constant", "local_linear"),
                    kernel=c("gaussian", "uniform"),
                    bandwidth=NULL, nmulti=5){
  trend <- match.arg(trend)
  kernel <- match.arg(kernel)
  x <- na.fail(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric")

  rmn <- rowMeans(x)
  if (trend == "grand_mean"){
    x <- x - mean(rmn)
  } else if (trend == "ensemble_means"){
    x <- x - rmn
  } else if (trend == "local_constant"
             || trend == "local_linear"){
    samplet <- as.integer(nrow(x))
    step <- seq(1, samplet)
    data <- data.frame(step=step, rmn=rmn)
    srm <- smooth(data=data, bandwidth=bandwidth, est=trend, kernel=kernel,
                  nmulti=nmulti)
    x <- x - srm$smooth
    bandwidth <- srm$bandwidth
  }
  list(x=x, bandwidth=bandwidth)
}

smooth <- function(data, est, bandwidth, kernel="gaussian", nmulti=5){
  is.constant <- grepl("constant", est)
  rt <- ifelse(is.constant, "lc", "ll")
  if (!is.null(bandwidth)){
    bw <- np::npregbw(formula=rmn ~ step, bws=bandwidth,
                      regtype=rt, ckertype=kernel,
                      bandwidth.compute=FALSE, data=data,
                      na.action=na.fail)
  } else {
    bw <- np::npregbw(formula=rmn ~ step, regtype=rt, ckertype=kernel,
                      bwmethod="cv.ls", bwtype="fixed",
                      data=data, na.action=na.fail, nmulti=nmulti)
  }
  mod <- np::npreg(bw)
  list(smooth=fitted(mod), bandwidth=bw$bw)
}

autocor <- function(x, cortype=c("correlation", "covariance"), lag=1,
                    bandwidth=NULL, trend=c("local_constant", "local_linear"),
                    kernel=c("gaussian", "uniform"), nmulti=5){
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
  smth_xx_lag <- smooth(data=data, bandwidth=bandwidth, kernel=kernel,
                est=trend, nmulti=nmulti)
  if (cortype == "correlation") {
    xx <- rowMeans(x1 * x1) * 0.5 + rowMeans(x2 * x2) * 0.5
    data <- data.frame(step=step, rmn=xx)
    smth_xx <- smooth(data=data, bandwidth=smth_xx_lag$bandwidth,
                      kernel=kernel, est=trend, nmulti=nmulti)
    list(smooth=smth_xx_lag$smooth / smth_xx$smooth,
         bandwidth=smth_xx$bandwidth)
  } else {
    smth_xx_lag
  }
}
