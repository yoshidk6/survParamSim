
#' @rdname calculate_hr_pi
#' @export
#'
calc_ave_hr_pi <- function(sim, trt, group = NULL, pi.range = 0.95,
                           calc.obs = TRUE, trt.assign = c("default", "reverse"),
                           time.max = 1000){

  # Replace nest with packageVersion("tidyr") >= '1.0.0' for a speed issue
  # and different behavior when no grouping is supplied
  # See https://github.com/tidyverse/tidyr/issues/751
  nest2 <- ifelse(utils::packageVersion("tidyr") == '1.0.0', tidyr::nest_legacy, tidyr::nest)
  unnest2 <- ifelse(utils::packageVersion("tidyr") == '1.0.0', tidyr::unnest_legacy, tidyr::unnest)

  trt.assign <- match.arg(trt.assign)

  # Handle trt variable -------------------------------------------------------------------

  handle.trt.group.output <- handle_trt_group(trt, trt.assign, group,
                                              sim$newdata.nona.obs, sim$newdata.nona.sim)
  group.syms       <- handle.trt.group.output$group.syms
  trt.sym          <- handle.trt.group.output$trt.sym
  newdata.nona.obs <- handle.trt.group.output$newdata.nona.obs
  newdata.nona.sim <- handle.trt.group.output$newdata.nona.sim
  trt.levels       <- handle.trt.group.output$trt.levels


  # Calc HR for observed data -------------------------------------------------------------------

  if(methods::is(sim, "survparamsim_pre_resampled")){
    if(sim$newdata.orig.missing & calc.obs) {
      warning("Original observed data not provided in `surv_param_sim_pre_resampled()` and HR will not be estimated for the observed data. Speficy `calc.obs = FALSE` to suppress this warning.")
      calc.obs = FALSE
    }
  }
  obs.hr <- calc_hr_for_obs(sim, newdata.nona.obs, group.syms, trt, trt.sym, trt.assign, trt.levels,
                            nest2, unnest2, calc.obs)


  # Calc HR for simulated data ----------------------------------------------------------------

  ## First nest data - cox fit will done for each nested data
  newdata.trt.group <-
    newdata.nona.sim %>%
    dplyr::select(subj.sim, !!trt.sym, !!!group.syms)

  sim.nested <-
    sim$sim %>%
    dplyr::left_join(newdata.trt.group, by = "subj.sim") %>%
    dplyr::group_by(rep, !!!group.syms) %>%
    nest2()


  # Calculate average HR
  ## sim.hr has rep, !!trt, !!!group, HR, description (= "sim") columns,
  ## This is what we need to emulate with the new method
  ## it also has p.value.coef.wald & p.value.logrank but it should be ok to ignore for now
  sim.hr <- calc_hr_with_average_surv(sim.nested, trt, trt.sym, trt.levels,
                                      time.max, unnest2)


  ## Reverse back the factor
  if(trt.assign == "reverse"){
    sim.hr <-
      sim.hr %>%
      dplyr::mutate(!!trt.sym := forcats::fct_rev(!!trt.sym))
  }


  # Calc quantiles ----------------------------------------------------------------

  ## Better to first implement 3 or more treatment groups handling
  ## before making the quantile function to work
  # hr.pi.quantile <- calc_hr_quantiles(pi.range, sim.hr, obs.hr, calc.obs,
  #                                     group.syms, trt.sym)

  # Output ---------------------------------------------------------------
  out <- list()

  out$calc.obs <- calc.obs
  out$pi.range   <- pi.range

  out$group.syms <- group.syms
  out$trt.sym    <- trt.sym

  out$obs.hr <- obs.hr
  out$sim.hr <- sim.hr
  # out$hr.pi.quantile <- hr.pi.quantile

  out$trt.levels <- trt.levels

  structure(out, class = c("survparamsim.hrpi"))
}



calc_hr_with_average_surv <- function(sim.nested, trt, trt.sym, trt.levels,
                                      time.max = 1000, unnest2) {

  ## Define function to calc HR
  calc_hr_each_sim <- function(x, y){

    # First level in the factor serves as the reference group
    lp.vec.control <- x$lp[as.numeric(x[[trt]]) == 1]
    lp.vec.treatment <- x$lp[as.numeric(x[[trt]]) == 2]

    ahr <-
      calc_ave_hr_from_lp(lp.vec.control, lp.vec.treatment, scale = y,
                          dist = "lognormal", time.max = time.max)
  }

  ## Calc HR
  sim.hr <-
    sim.nested %>%
    dplyr::ungroup() %>%
    dplyr::mutate(HR = purrr::map2_dbl(data, sim$scale.vec, calc_hr_each_sim),
                  description = "sim") %>%
    dplyr::select(-data)

  return(sim.hr)
}



#' Calculate average hazard ratio from linear predictors
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
calc_ave_hr_from_lp <- function(lp.vec.control, lp.vec.treatment, scale = NULL,
                                dist = "lognormal",
                                time.max = 1000){

  # Currently only accept log normal
  dist <- match.arg(dist)

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
create_survfun <- function(lpvec, scale, dist = "lognormal"){
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
create_pdf <- function(lpvec, scale, dist = "lognormal"){
  function(x) {
    pdmatrix <-
      vapply(lpvec,
             function(lp) {dlnorm(x=x, meanlog=lp, sdlog=exp(scale))},
             numeric(length(x)))

    if(length(x) == 1) return(mean(pdmatrix))
    return(rowMeans(pdmatrix))
  }
}


