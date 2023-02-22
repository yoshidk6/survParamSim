
#' Calculate average hazard ratio
#'
#' Survival and PDF functions used in calculation are average of these functions for subjects in
#' the individual groups, because every subject has different survival and PDF functions
#' @param lp.vec.control A vector of linear predictor (lp) for the control (or reference) group.
#' The number of elements equal to the number of subjects in the control group of `newdata`
#'  used for `surv_param_sim()` function.
#' @param lp.vec.treatment A vector of linear predictor (lp) for the treatment (or test) group.
#' @param scale Scale variable used for simulation (NULL for exponential model)
#' @dist Distribution for the parametric survival model
#' @time.max Time for calculation of average HR
#' @return Average hazard ratio from time `0` to `time.max`
calc_ave_hr_oc <- function(lp.vec.control, lp.vec.treatment, scale = NULL,
                           dist = "lognormal",
                           time.max = 1000){

  survfun.control   <- create_survfun(lp.vec.control,   scale, dist = dist)
  survfun.treatment <- create_survfun(lp.vec.treatment, scale, dist = dist)
  pdf.control   <- create_pdf(lp.vec.control,   scale, dist = dist)
  pdf.treatment <- create_pdf(lp.vec.treatment, scale, dist = dist)

  integrand1 <- function(x){survfun.control(x) * pdf.treatment(x)}
  integrand2 <- function(x){survfun.treatment(x) * pdf.control(x)}

  term1 <- integrate(integrand1, lower=0, upper=time.max)$value
  term2 <- integrate(integrand2, lower=0, upper=time.max)$value

  ahr <- term1 / term2

  return(ahr)
}




#' Create a function to calculate S(t)
#' The function to be created will take time (`x`) as the input argument and return
#'
create_survfun <- function(lpvec, scale){
  function(x) {
    survmatrix <-
      vapply(lpvec,
             function(lp) {plnorm(q=x, meanlog=lp, sdlog=exp(scale), lower=FALSE)},
             numeric(length(x)))

    if(length(x) == 1) return(mean(survmatrix))
    return(rowMeans(survmatrix))
  }
}

# Create a function to calculate PDF(t)
create_pdf <- function(lpvec, scale){
  function(x) {
    pdmatrix <-
      vapply(lpvec,
             function(lp) {dlnorm(x=x, meanlog=lp, sdlog=exp(scale))},
             numeric(length(x)))

    if(length(x) == 1) return(mean(pdmatrix))
    return(rowMeans(pdmatrix))
  }
}
