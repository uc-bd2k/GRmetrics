
.onAttach <- function(...) {
    info <- c(
      "Try our web application at http://www.grcalculator.org",
      "[other messages here]"
    )
    packageStartupMessage(paste(strwrap(info), collapse = "\n"))
}