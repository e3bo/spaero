
#' Detrend ensemble time series
#'
#' @param x A univariate or multivariate numeric time series object or
#'   a numeric vector or matrix.
#' @param lag Lag at which to calculate the acf. Default is 1. This
#'   lag is in terms of the index of \code{x} and does not account for the
#'   frequency of \code{x} if \code{x} is a time series.
#' @param type Character string giving the type of acf to be
#'   returned. Allowed values are '"correlation"' (the default), and
#'   '"covariance"'. Will be partially matched.
#' @param type Character string giving method of calculating the trend
#' to subtract. Allowed values are '"grand.mean"' (the default),
#' '"ensemble.means"', '"local.mean.gaussian.weights"',
#' '"local.linear.gaussian.weights"', '"none"'. Will be partially
#' matched.
#' @param bandwidth Bandwith for any kernel-based detrending
#'   done. If not supplied, a bandwidth will be selected by
#'   cross-validation of least squared error.
#'
#' @details Any missing values in 'x' will cause an error.
#' @export
#'
detrend <- function(x, type=c("grand.mean", "ensemble.means",
                        "local.constant.gaussian.weights",
                        "local.linear.gaussian.weights"),
                    bandwidth=NULL, nmulti=5){
  type <- match.arg(type)
  x <- na.fail(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric")

  rmn <- rowMeans(x)
  if (type == "grand.mean"){
    x <- x - mean(rmn)
  } else if (type == "ensemble.means"){
    x <- x - rmn
  } else if (type == "local.constant.gaussian.weights"
             || type == "local.linear.gaussian.weights"){
    samplet <- as.integer(nrow(x))
    step <- seq(1, samplet)
    data <- data.frame(step=step, rmn=rmn)
    srm <- smooth(data, type, bandwidth, nmulti)
    x <- x - srm
  }
  x
}

smooth <- function(data, type, bandwidth, nmulti=5){
  is.constant <- grepl("constant", type)
  rt <- ifelse(is.constant, "lc", "ll")
  if (!is.null(bandwidth)){
    bw <- np::npregbw(formula=rmn ~ step, bws=bandwidth,
                      regtype=rt, ckertype="gaussian", ckerorder=2,
                      bandwidth.compute=FALSE, data=data,
                      na.action=na.fail)
  } else {
    bw <- np::npregbw(formula=rmn ~ step, regtype=rt, ckertype="gaussian",
                      bkerorder=2, bwmethod="cv.ls", bwtype="fixed",
                      data=data, na.action=na.fail, nmulti=nmulti)
  }
  mod <- np::npreg(bw)
  fitted(mod)
}

autocor <- function(x, coretype=c('correlation', 'covariance'), lag=1,
                    bandwidth=NULL,
                    smoothtype=c("local.constant.gaussian.weights",
                        "local.linear.gaussian.weights"), nmulti=5){
  smoothtype <- match.arg(smoothtype)
  coretype <- match.arg(coretype)
  x <- na.fail(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric")
  stopifnot(lag >= 0)

  n <- nrow(x)
  end1 <- n - lag
  start2 <- 1 + lag
  x1 <- x[1:end1, , drop=FALSE]
  x2 <- x[start2:n, , drop=FALSE]
  xx <- rowMeans(x1 * x2)

  step <- seq(start2, n)
  data <- data.frame(step=step, rmn=xx)
  sxx <- smooth(data, smoothtype, bandwidth, nmulti)
  sxx
}
