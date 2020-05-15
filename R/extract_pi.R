#' Functions to extract prediction intervals and observed data
#'
#'
#' @rdname extractpi
#' @export
#' @param km.pi A return object from \code{\link{calc_km_pi}} function.
#' @param trunc.sim.censor A logical specifying whether to truncate the simulated
#' curve at the last time of `censor.dur`` specified in \code{\link{surv_param_sim}}.
#' @details
#' \code{\link{extract_km_pi}} extracts prediction intervals of simulated Kaplan-Meier curves.
extract_km_pi <- function(km.pi, trunc.sim.censor = TRUE) {

  group <- km.pi$group
  trt   <- km.pi$trt

  obs.km <- km.pi$obs.km
  sim.km.quantile <- km.pi$sim.km.quantile

  group.syms <- rlang::syms(group)
  trt.syms   <- rlang::syms(trt)


  #### Below will not work if simtimelast is missing and obs.km is not calculated ####


  # Extract quantile from simulation
  # Limit data based on `simtimelast` or the last observed time
  if(is.null(km.pi$simtimelast) & km.pi$calc.obs){
    ## Get last obs time for each group
    timelast <-
      obs.km %>%
      dplyr::group_by(!!!trt.syms, !!!group.syms) %>%
      dplyr::arrange(time) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!!trt.syms, !!!group.syms, timelast = time)

    if(is.null(c(group, trt))){
      sim.km.quantile.plot <-
        sim.km.quantile %>%
        tidyr::crossing(timelast) %>%
        dplyr::filter(time <= timelast) %>%
        dplyr::select(-timelast)

    } else {
      sim.km.quantile.plot <-
        sim.km.quantile %>%
        dplyr::full_join(timelast, by = c(group, trt)) %>%
        dplyr::filter(time <= timelast) %>%
        dplyr::select(-timelast)
    }

  } else if(!is.null(km.pi$simtimelast)) {
    sim.km.quantile.plot <-
      sim.km.quantile %>%
      dplyr::filter(time <= km.pi$simtimelast)

  } else {
    sim.km.quantile.plot <-
      sim.km.quantile %>%
      dplyr::filter(time <= km.pi$t.last)
  }


  if(trunc.sim.censor & !is.null(km.pi$censor.dur)){
    sim.km.quantile.plot <-
      sim.km.quantile.plot %>%
      dplyr::filter(time <= km.pi$censor.dur[[2]])

  }

  return(sim.km.quantile.plot)
}


#' @rdname extractpi
#' @export
#' @param outtype Specifies whether output will be in long or wide format.
#' @details
#' \code{\link{extract_median_surv}} extracts prediction intervals of
#' median survival times and and the corresponding observed values.
extract_median_surv <- function(km.pi, outtype = c("long", "wide")) {

  outtype <- match.arg(outtype)

  out <- km.pi$median.pi

  if(outtype == "wide"){
    out <-
      out %>%
      dplyr::select(-quantile) %>%
      tidyr::spread(description, median)

    if(km.pi$calc.obs){
      out <- dplyr::select(out, pi_low, pi_med, pi_high, obs, dplyr::everything())
    } else {
      out <- dplyr::select(out, pi_low, pi_med, pi_high, dplyr::everything())
    }
  }

  return(out)
}

#' @rdname extractpi
#' @export
#' @param hr.pi a return object from \code{\link{calc_hr_pi}} function.
#' @details
#' \code{\link{extract_hr_pi}} extracts prediction intervals of simulated
#' hazard ratios and the corresponding observed values.
extract_hr_pi <- function(hr.pi, outtype = c("long", "wide")) {

  outtype <- match.arg(outtype)

  out <- hr.pi$hr.pi.quantile

  if(outtype == "wide"){
    out <-
      out %>%
      dplyr::select(-quantile) %>%
      tidyr::spread(description, HR)

    if(hr.pi$calc.obs){
      out <- dplyr::select(out, pi_low, pi_med, pi_high, obs, dplyr::everything())
    } else {
      out <- dplyr::select(out, pi_low, pi_med, pi_high, dplyr::everything())
    }
  }

  return(out)

}

