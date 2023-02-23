
#' @rdname calculate_hr_pi
#' @export
#'
#' @param simtimelast An optional numeric to specify duration for average HR calculation.
#' If NULL (default), the last observation time in the `newdata` will be used.
#' @param boot.subj Boolean to specify whether bootstrapping of subjects are performed
#' before calculating HR. Default TRUE.
calc_ave_hr_pi <- function(sim, trt, group = NULL, pi.range = 0.95,
                           calc.obs = TRUE, trt.assign = c("default", "reverse"),
                           simtimelast = NULL, boot.subj = TRUE){

  trt.assign <- match.arg(trt.assign)


  ## time for output
  if(is.null(simtimelast)) simtimelast <- sim$t.last.orig.new

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
                            calc.obs)


  # Calc HR for simulated data ----------------------------------------------------------------

  ## First nest data - cox fit will done for each nested data
  newdata.trt.group <-
    newdata.nona.sim %>%
    dplyr::select(subj.sim, !!trt.sym, !!!group.syms)

  sim.nested.2 <-
    sim$sim %>%
    dplyr::left_join(newdata.trt.group, by = "subj.sim") %>%
    dplyr::group_by(rep, !!trt.sym, !!!group.syms) %>%
    tidyr::nest() %>%
    dplyr::ungroup()

  if (boot.subj) {
    sim.nested.2 <-
      sim.nested.2 %>%
      dplyr::mutate(data = purrr::map(data, function(x) dplyr::slice_sample(x, prop = 1, replace = TRUE)))
  }

  # Calculate average HR
  sim.hr <- calc_hr_with_average_surv(sim.nested.2, sim$scale.bs.df,
                                      trt, trt.sym, trt.levels, group, group.syms,
                                      dist = sim$survreg$dist,
                                      simtimelast)

  ## Reverse back the factor
  if(trt.assign == "reverse"){
    sim.hr <-
      sim.hr %>%
      dplyr::mutate(!!trt.sym := forcats::fct_rev(!!trt.sym))
  }


  # Calc quantiles ----------------------------------------------------------------

  # Better to first implement 3 or more treatment groups handling
  # before making the quantile function to work
  hr.pi.quantile <- calc_hr_quantiles(pi.range, sim.hr, obs.hr, calc.obs,
                                      group.syms, trt.sym)

  # Output ---------------------------------------------------------------
  out <- list()

  out$calc.obs <- calc.obs
  out$pi.range   <- pi.range

  out$group.syms <- group.syms
  out$trt.sym    <- trt.sym

  out$obs.hr <- obs.hr
  out$sim.hr <- sim.hr
  out$hr.pi.quantile <- hr.pi.quantile

  out$trt.levels <- trt.levels

  structure(out, class = c("survparamsim.hrpi.aveHR", "survparamsim.hrpi"))
}



calc_hr_with_average_surv <- function(sim.nested.2, scale.bs.df,
                                      trt, trt.sym, trt.levels, group, group.syms, dist,
                                      simtimelast = 1000) {

  # Extract linear predictor (lp), also get scale
  df.lp.extracted <-
    sim.nested.2 %>%
    dplyr::mutate(lp = purrr::map(data, function(x) x$lp)) %>%
    dplyr::select(-data) %>%
    dplyr::left_join(scale.bs.df, by = "rep")

  df.lp.control <-
    df.lp.extracted %>%
    dplyr::filter(as.numeric(!!trt.sym) == 1) %>%
    dplyr::select(-!!trt.sym)
  df.lp.treatment <-
    df.lp.extracted %>%
    dplyr::filter(as.numeric(!!trt.sym) != 1)

  # Create survival and PDF functions
  df.surv.pdf.fun.control <-
    df.lp.control %>%
    dplyr::mutate(survfun.control = purrr::map2(lp, scale, function(x, y) create_survfun(lpvec = x, scale = y, dist = dist)),
                  pdf.control     = purrr::map2(lp, scale, function(x, y) create_pdf(lpvec = x, scale = y, dist = dist))) %>%
    dplyr::select(-lp, -scale)
  df.surv.pdf.fun.treatment <-
    df.lp.treatment %>%
    dplyr::mutate(survfun.trt = purrr::map2(lp, scale, function(x, y) create_survfun(lpvec = x, scale = y, dist = dist)),
                  pdf.trt     = purrr::map2(lp, scale, function(x, y) create_pdf(lpvec = x, scale = y, dist = dist))) %>%
    dplyr::select(-lp, -scale)

  df.surv.pdf.fun.join <-
    dplyr::left_join(df.surv.pdf.fun.treatment,
                     df.surv.pdf.fun.control,
                     by = c("rep", group))

  # Calculate aveHR
  sim.hr <-
    df.surv.pdf.fun.join %>%
    dplyr::mutate(integrand1 = purrr::map2(survfun.control, pdf.trt, function(x, y) function(t){x(t) * y(t)}),
                  integrand2 = purrr::map2(survfun.trt, pdf.control, function(x, y) function(t){x(t) * y(t)})) %>%
    dplyr::mutate(term1 = purrr::map_dbl(integrand1, function(x) integrate(x, lower = 0, upper = simtimelast)$value),
                  term2 = purrr::map_dbl(integrand2, function(x) integrate(x, lower = 0, upper = simtimelast)$value)) %>%
    dplyr::mutate(HR = term1 / term2) %>%
    dplyr::select(rep, !!!group.syms, !!trt.sym, HR) %>%
    dplyr::mutate(description = "sim")

  return(sim.hr)
}



# Calculate average hazard ratio from linear predictors
#
# Only used for code testing only, as this workflow was not the best fit for
# calculating HR with 3 or more treatment arms
#
# Survival and PDF functions used in calculation are average of these functions for subjects in
# the individual groups, because every subject has different survival and PDF functions
calc_ave_hr_from_lp <- function(lp.vec.control, lp.vec.treatment, scale = NULL,
                                dist = "lognormal",
                                simtimelast = NULL){

  # Currently only accept log normal
  dist <- match.arg(dist)

  survfun.control   <- create_survfun(lp.vec.control,   scale, dist = dist)
  survfun.treatment <- create_survfun(lp.vec.treatment, scale, dist = dist)
  pdf.control   <- create_pdf(lp.vec.control,   scale, dist = dist)
  pdf.treatment <- create_pdf(lp.vec.treatment, scale, dist = dist)

  integrand1 <- function(x){survfun.control(x) * pdf.treatment(x)}
  integrand2 <- function(x){survfun.treatment(x) * pdf.control(x)}

  term1 <- integrate(integrand1, lower=0, upper=simtimelast)$value
  term2 <- integrate(integrand2, lower=0, upper=simtimelast)$value

  ahr <- term1 / term2

  return(ahr)
}



# Create a function to calculate S(t)
# The function to be created will take time (`x`) as the input argument and return
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


