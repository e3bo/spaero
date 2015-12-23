if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free()
  })
}

set.seed(123)

context("detrending")

test_that("Mean-based detrending works", {
  expect_equal(detrend(1:10, trend="grand.mean")$x, matrix(1:10 - 5.5))
  expect_equal(detrend(cbind(1:10, 2:11), trend="grand.mean")$x,
               cbind(1:10, 2:11) - 6)
  expect_equal(detrend(1:10, trend="ensemble")$x, matrix(rep(0, 10)))
  expect_equal(detrend(cbind(1:10, 2:11), trend="ensemble")$x,
               cbind(1:10, 2:11) - 1:10 - 0.5)
  })

context("smoothing")

test_that("Smoothing function works as expected", {
  simple_smooth <- function(data, est, kernel="gaussian", bandwidth){
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
    if (est == "local.constant"){
      colSums(w * data$rmn)
    } else {
      wlm <- function(x){
        m <- lm(rmn~step, weights=w[, x], data=data)
        fitted(m)[x]
      }
      sapply(seq_along(data$step), wlm)
    }
  }
  data <- data.frame(step=1:10, rmn=1:10)
  expect_equal(smooth(data, est="local.constant", bandwidth=2)$smooth,
               simple_smooth(data, est="local.constant",  bandwidth=2))
  expect_equal(smooth(data, est="local.linear", bandwidth=2)$smooth,
               simple_smooth(data, est="local.linear",  bandwidth=2),
               check.names=FALSE)

  expect_equal(smooth(data, est="local.constant", kernel="uniform",
                      bandwidth=3)$smooth,
               simple_smooth(data, est="local.constant", kernel="uniform",
                             bandwidth=3))
  expect_equal(smooth(data, est="local.linear", bandwidth=3)$smooth,
               simple_smooth(data, est="local.linear", bandwidth=3),
               check.names=FALSE)

  data2 <- data.frame(step=20:100, rmn=rnorm(81))
  expect_equal(smooth(data2, est="local.constant", bandwidth=5)$smooth,
               simple_smooth(data2, est="local.constant", bandwidth=5))
  expect_equal(smooth(data2, est="local.linear", bandwidth=5,
                      kernel="uniform")$smooth,
               simple_smooth(data2, est="local.linear", bandwidth=5,
                             kernel="uniform"),
               check.names=FALSE)
})
