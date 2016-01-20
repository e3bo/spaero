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
    if (est == "local_constant"){
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
  expect_equal(smooth(data, est="local_constant", bandwidth=2)$smooth,
               simple_smooth(data, est="local_constant",  bandwidth=2))
  expect_equal(smooth(data, est="local_linear", bandwidth=2)$smooth,
               simple_smooth(data, est="local_linear",  bandwidth=2),
               check.names=FALSE)

  expect_equal(smooth(data, est="local_constant", kernel="uniform",
                      bandwidth=3)$smooth,
               simple_smooth(data, est="local_constant", kernel="uniform",
                             bandwidth=3))
  expect_equal(smooth(data, est="local_linear", bandwidth=3)$smooth,
               simple_smooth(data, est="local_linear", bandwidth=3),
               check.names=FALSE)

  data2 <- data.frame(step=20:100, rmn=rnorm(81))
  expect_equal(smooth(data2, est="local_constant", bandwidth=5)$smooth,
               simple_smooth(data2, est="local_constant", bandwidth=5))
  expect_equal(smooth(data2, est="local_linear", bandwidth=5,
                      kernel="uniform")$smooth,
               simple_smooth(data2, est="local_linear", bandwidth=5,
                             kernel="uniform"),
               check.names=FALSE)
})

context("cross validation")

test_that("Cross validation of mean bandwidth produces reasonable results", {
     w <- rnorm(100)
     xnext <- function(xlast, w) 0.1 * xlast + w
     x <- Reduce(xnext, x=w, init=0, accumulate=TRUE)
     dacf <- get_dynamic_acf(x, center_trend="local_constant",
                             acf_bandwidth=length(x))
     expect_gt(dacf$centered$bandwidth, length(x) / 2)
     xx <- x + seq_along(x) * 100
     dacf <- get_dynamic_acf(xx, center_trend="local_constant",
                             acf_bandwidth=length(x))
     expect_lt(dacf$centered$bandwidth, 1)
     xx <- x + sin(seq_along(x) / 10) * 10
     dacf <- get_dynamic_acf(xx, center_trend="local_linear",
                             acf_bandwidth=length(x))
     expect_gt(dacf$centered$bandwidth, 1)
     expect_lt(dacf$centered$bandwidth, 5)
})

test_that("Cross validation of acf bandwidth produces reasonable results", {
     n <- 100
     gen_ar1 <- function(n, phi=0.1){
       w <- rnorm(n)
       init <- rnorm(n=1, sd=sqrt(1 / (1 - phi^2)))
       xnext <- function(xlast, w) phi * xlast + w
       Reduce(xnext, x=w, init=init, accumulate=TRUE)
     }
     set.seed(123)
     n <- 100
     nts <- 100
     x <- replicate(nts, gen_ar1(n=n, phi=0.1))
     dacf <- get_dynamic_acf(x, center_trend="local_linear",
                             center_bandwidth=n,
                             acf_trend="local_linear")
     expect_gt(dacf$acf$bandwidth, n / 2)

     nts <- 100
     n <- 100
     x <- replicate(nts, gen_ar1(n=n, phi=0.9))
     dacf <- get_dynamic_acf(x, center_trend="local_linear",
                             center_bandwidth=n,
                             acf_trend="local_linear")

})
