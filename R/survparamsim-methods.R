#' Methods for S3 objects in the package
#'
#' @name survparamsim-methods
#'
NULL


#' @rdname survparamsim-methods
#' @export
#' @param x An object of the corresponding class
#' @param object An object of the corresponding class
#' @param ... Additional arguments passed to methods.
print.survparamsim <- function(x, ...) {
  cat("---- Simulated survival data with the following model ----\n")
  dput(x$survreg$call)
  cat("\n")
  cat("* Use `extract_sim()` function to extract individual simulated survivals\n")
  cat("* Use `calc_km_pi()` function to get survival curves and median survival time\n")
  cat("* Use `calc_hr_pi()` function to get hazard ratio\n\n")
  cat("* Settings:\n")
  cat("    #simulations:", max(x$sim$rep), "\n", sep=" ")
  cat("    #subjects:", x$newdata.nona.obs %>% nrow(),
      "(without NA in model variables)\n", sep=" ")
}

