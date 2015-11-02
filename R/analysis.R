
#' Estimate the autocorrelation at a given lag over rolling windows
#'
#' @param x A univariate or multivariate numeric time series object or
#'   a numeric vector or matrix.
#' @param lag Lag at which to calculate the acf. Default is 1. This
#'   lag is in terms of the index of \code{x} and does not account for the
#'   frequency of \code{x} if \code{x} is a time series.
#' @param type Character string giving the type of acf to be
#'   returned. Allowed values are '"correlation"' (the default), and
#'   '"covariance"'. Will be partially matched.
#' @param na.action Function to be called to handle missing
#'   values. 'na.pass' can be used.
#' @param detrend Character string giving method of centering the
#'   series before calculating the ACF. Allowed values are '"grand.mean"'
#'   (the default), '"ensemble.means"', '"local.mean.gaussian.weights"',
#'   '"local.linear.gaussian.weights"', '"none"'. Will be partially matched.
#' @param bandwidth.detrend Bandwith for any kernal-based detrending
#'   done. If not supplied, a bandwidth will be selected by
#'   cross-validation of least squared error.
#' @export
#'
#'
GetRollingACF <- function(x, lag=1, type=c("correlation", "covariance"),
                          na.action=na.fail, detrend=c("grand.mean",
                                                 "ensemble.means",
                                                 "local.constant.gaussian.weights",
                                                 "local.linear.gaussian.weights",
                                                 "none"),
                          bandwidth.detrend=NULL){
  ## First few lines from stats::acf() function
  detrend <- match.arg(detrend)
  type <- match.arg(type)
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric")
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  rmn <- rowMeans(x, na.rm=TRUE)
  step <- seq(1, sampleT)
  if (detrend == "grand.mean"){
    x <- x - mean(rmn)
  } else if (detrend == "ensemble.means"){
    x <- x - rmn
  } else if (detrend == "local.constant.gaussian.weights"
             || detrend == "local.linear.gaussian.weights"){
    data <- data.frame(step=step, rmn=rmn)
    is.constant <- grepl('constant', detrend)
    rt <- ifelse(is.constant, 'lc', 'll')
    if (!is.null(bandwidth.detrend)){
      bw <- np::npregbw(formula=rmn ~ step, bws=bandwidth.detrend,
                        regtype=rt, ckertype="gaussian", ckerorder=2,
                        bandwidth.compute=FALSE, data=data)
    } else {
      bw <- np::npregbw(formula=rmn ~ step, regtype=rt, ckertype="gaussian",
                        bkerorder=2, bwmethod='cv.ls', bwtype='fixed',
                        data=data)
    }
    mod <- np::npreg(bw)
    srm <- fitted(mod)
    x <- x - srm
  }
  x
}
