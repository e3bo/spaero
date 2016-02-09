if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free()
  })
}

set.seed(123)

context("detrending")

test_that("Mean-based detrending works", {
  expect_equal(detrend(1:10, trend="grand_mean")$x, matrix(1:10 - 5.5))
  expect_equal(detrend(cbind(1:10, 2:11), trend="grand_mean")$x,
               cbind(1:10, 2:11) - 6)
  expect_equal(detrend(1:10, trend="ensemble")$x, matrix(rep(0, 10)))
  expect_equal(detrend(cbind(1:10, 2:11), trend="ensemble")$x,
               cbind(1:10, 2:11) - 1:10 - 0.5)
})

test_that("Skipping detrending works", {
  expect_equal(detrend(1:10, trend="assume_zero")$x, matrix(1:10))
})

if (requireNamespace("np", quietly = TRUE)) {
  if (is.null(options("np.messages")$np.messages)) options(np.messages = TRUE)
  if (is.null(options("np.tree")$np.tree)) options(np.tree = FALSE)

  context("smoothing")

  test_that("Smoothing function works as expected", {
    np_smooth <- function(data, est, bandwidth, kernel="gaussian"){
      is.constant <- grepl("constant", est)
      rt <- ifelse(is.constant, "lc", "ll")
      bw <- np::npregbw(formula=rmn ~ step, bws=bandwidth,
                        regtype=rt, ckertype=kernel,
                        bandwidth.compute=FALSE, data=data,
                        na.action=na.fail)
      mod <- np::npreg(bw)
      list(smooth=fitted(mod), bandwidth=bw$bw)
    }
    data <- data.frame(step=1:10, rmn=1:10)
    expect_equal(smooth(data, est="local_constant", bandwidth=2)$smooth,
                 np_smooth(data, est="local_constant",  bandwidth=2)$smooth)
    expect_equal(smooth(data, est="local_linear", bandwidth=2)$smooth,
                 np_smooth(data, est="local_linear",  bandwidth=2)$smooth,
                 check.names=FALSE)

    expect_equal(smooth(data, est="local_constant", kernel="uniform",
                        bandwidth=3)$smooth,
                 np_smooth(data, est="local_constant", kernel="uniform",
                           bandwidth=3)$smooth)
    expect_equal(smooth(data, est="local_linear", bandwidth=3)$smooth,
                 np_smooth(data, est="local_linear", bandwidth=3)$smooth,
                 check.names=FALSE)

    data2 <- data.frame(step=20:100, rmn=rnorm(81))
    expect_equal(smooth(data2, est="local_constant", bandwidth=5)$smooth,
                 np_smooth(data2, est="local_constant", bandwidth=5)$smooth)
    expect_equal(smooth(data2, est="local_linear", bandwidth=5,
                        kernel="uniform")$smooth,
                 np_smooth(data2, est="local_linear", bandwidth=5,
                           kernel="uniform")$smooth,
                 check.names=FALSE)
    })
}

context("smoothing arguments")

test_that("invalid bandwidths lead to errors", {
  data <- data.frame(step=1:10, rmn=1:10)
  expect_error(smooth(data),
    regexp="argument \"bandwidth\" is missing, with no default")
  expect_error(smooth(data, bandwidth=c(1:20)),
    regexp="argument \"bandwidth\" must be provided as a single numeric value")
  expect_error(smooth(data, bandwidth="hmm"),
    regexp="argument \"bandwidth\" must be provided as a single numeric value")
  expect_error(smooth(data, bandwidth=0.5),
    regexp="argument \"bandwidth\" must be >= 1")
          })
