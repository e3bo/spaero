.onLoad <- function (lib, pkg) {
  if(is.null(options("np.messages")$np.messages))
    options(np.messages = TRUE)
  if(is.null(options("np.tree")$np.tree))
    options(np.tree = FALSE)
}
