if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free()
  })
}

set.seed(123)

context("detrending")

test_that("Mean-based detrending works", {
  expect_equal(detrend(1:10, type="grand.mean"), matrix(1:10 - 5.5))
  expect_equal(detrend(cbind(1:10, 2:11), type="grand.mean"),
               cbind(1:10, 2:11) - 6)
  expect_equal(detrend(1:10, type="ensemble"), matrix(rep(0, 10)))
  expect_equal(detrend(cbind(1:10, 2:11), type="ensemble"),
               cbind(1:10, 2:11) - 1:10 - 0.5)
  })

context("smoothing")

test_that("Smoothing function works as expected", {
  simple_smooth <- function(data, type, bandwidth){
    kern <- function(ind, bw=bandwidth){
      dist <- abs(data$step - ind) / bw
      w <- dnorm(dist)
      w / sum(w)
    }
    w <- sapply(data$step, kern)
    if (type == "local.constant"){
      colSums(w * data$rm)
    } else {
      wlm <- function(x){
        m <- lm(rmn~step, weights=w[, x], data=data)
        fitted(m)[x]
      }
      sapply(seq_along(data$step), wlm)
    }
  }
  data <- data.frame(step=1:10, rmn=1:10)
  expect_equal(smooth(data, type="local.constant", bandwidth=2),
               simple_smooth(data, type="local.constant",  bandwidth=2))
  expect_equal(smooth(data, type="local.linear", bandwidth=2),
               simple_smooth(data, type="local.linear",  bandwidth=2),
               check.names=FALSE)

  expect_equal(smooth(data, type="local.constant", bandwidth=3),
               simple_smooth(data, type="local.constant", bandwidth=3))
  expect_equal(smooth(data, type="local.linear", bandwidth=3),
               simple_smooth(data, type="local.linear", bandwidth=3),
               check.names=FALSE)

  data2 <- data.frame(step=20:100, rmn=rnorm(81))
  expect_equal(smooth(data2, type="local.constant", bandwidth=5),
               simple_smooth(data2, type="local.constant", bandwidth=5))
  expect_equal(smooth(data2, type="local.linear", bandwidth=5),
               simple_smooth(data2, type="local.linear", bandwidth=5),
               check.names=FALSE)

})
