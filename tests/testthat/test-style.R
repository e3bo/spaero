context("lints")

test_that("Package Style", {
  skip_if_not_installed("lintr")
  skip_on_cran()
  lintr::expect_lint_free()
})
