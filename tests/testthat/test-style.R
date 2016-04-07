
has_lintr <- function() requireNamespace("lintr", quietly = TRUE)

context("lints")

test_that("Package Style", {
  if (!has_lintr()) skip("lintr not available")
  lintr::expect_lint_free()
})


