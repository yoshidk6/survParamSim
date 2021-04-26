#' Functions to extract prediction intervals and observed data
#'
#' @name extractpi
NULL


#' @rdname extractpi
#' @export
#' @param hr.pi a return object from [calc_hr_pi()] function.
#' @details
#' [extract_hr_pi()] extracts prediction intervals of simulated
#' hazard ratios and the corresponding observed values.
extract_hr_pi <- function(hr.pi, outtype = c("long", "wide")) {

  outtype <- match.arg(outtype)

  out <- hr.pi$hr.pi.quantile
  group.syms <- hr.pi$group.syms
  trt.sym    <- hr.pi$trt.sym

  if(outtype == "wide"){
    out <-
      out %>%
      dplyr::select(-quantile) %>%
      tidyr::spread(description, HR)

    if(hr.pi$calc.obs){
      out <- dplyr::select(out, !!!group.syms, !!trt.sym, pi_low, pi_med, pi_high, obs)
    } else {
      out <- dplyr::select(out, !!!group.syms, !!trt.sym, pi_low, pi_med, pi_high)
    }
  }

  return(out)

}


#' @rdname extractpi
#' @export
#' @param km.pi A return object from [calc_km_pi()] function.
#' @param trunc.sim.censor A logical specifying whether to truncate the simulated
#' curve at the last time of `censor.dur` specified in [surv_param_sim()].
#' @details
#' [extract_km_pi()] extracts prediction intervals of simulated Kaplan-Meier curves.
extract_km_pi <- function(km.pi, trunc.sim.censor = TRUE) {

  obs.km <- km.pi$obs.km
  sim.km.quantile <- km.pi$sim.km.quantile

  group.syms <- km.pi$group.syms
  trt.syms   <- km.pi$trt.syms


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

    if(length(c(group.syms, trt.syms)) == 0){
      sim.km.quantile.plot <-
        sim.km.quantile %>%
        tidyr::crossing(timelast) %>%
        dplyr::filter(time <= timelast) %>%
        dplyr::select(-timelast)

    } else {
      sim.km.quantile.plot <-
        sim.km.quantile %>%
        dplyr::full_join(timelast, by = as.character(c(group.syms, trt.syms))) %>%
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
#' [extract_medsurv_pi()] extracts prediction intervals of
#' median survival times and and the corresponding observed values.
extract_medsurv_pi <- function(km.pi, outtype = c("long", "wide")) {

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
#' @param outtype Specifies whether output will be in long or wide format.
#' @details
#' [extract_medsurv_delta_pi()] extracts prediction intervals of
#' delta of median survival times between treatment groups
extract_medsurv_delta_pi <- function(km.pi, outtype = c("long", "wide")) {

  pi.range   <- km.pi$pi.range
  group.syms <- km.pi$group.syms
  trt.sym    <- km.pi$trt.syms[[1]]

  outtype <- match.arg(outtype)

  sim.median.time.delta <- extract_medsurv_delta(km.pi)

  out <-
    sim.median.time.delta %>%
    dplyr::group_by(!!!group.syms, !!trt.sym) %>%
    dplyr::summarize(pi_low = as.numeric(stats::quantile(median_delta, probs = 0.5 - pi.range/2, na.rm = TRUE)),
                     pi_med = as.numeric(stats::quantile(median_delta, probs = 0.5, na.rm = TRUE)),
                     pi_high= as.numeric(stats::quantile(median_delta, probs = 0.5 + pi.range/2, na.rm = TRUE)),
                     .groups = "drop")

  quantiles <-
    tibble::tibble(description = c("pi_low", "pi_med", "pi_high"),
                   quantile = c(0.5 - pi.range/2, 0.5, 0.5 + pi.range/2))

  if(outtype == "long"){
    out <-
      out %>%
      tidyr::gather(description, median_delta, pi_low:pi_high) %>%
      dplyr::left_join(quantiles, by = "description") %>%
      dplyr::arrange(!!!group.syms, !!trt.sym, quantile)
  }

  return(out)
}





#' Functions to extract prediction intervals and observed data
#'
#' \lifecycle{deprecated}
#'
#' [extract_median_surv()] was renamed to [extract_medsurv_pi()] for function name consistency.
#'
#' @rdname extractpi_deprecated
#' @export
#' @param km.pi A return object from [calc_km_pi()] function.
#' @param outtype Specifies whether output will be in long or wide format.
#' @details
#' [extract_median_surv()] extracts prediction intervals of
#' median survival times and and the corresponding observed values.
extract_median_surv <- function(km.pi, outtype = c("long", "wide")) {
  lifecycle::deprecate_warn("0.1.5", "extract_median_surv()", "extract_medsurv_pi()")

  return(extract_medsurv_pi(km.pi, outtype))
}
