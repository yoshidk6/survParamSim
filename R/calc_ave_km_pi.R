
#' @rdname calculate_km_pi
#' @export
#'
#' @param calc.median.surv Whether to calculate median survival time for
#' [calc_ave_km_pi()]. Default FALSE as the calculation can be long.
#' Currently median survival calculation not implemented yet.
#'
calc_ave_km_pi <- function(sim, trt=NULL, group=NULL, pi.range = 0.95,
                           calc.obs = TRUE, simtimelast = NULL,
                           trt.assign = c("default", "reverse"),
                           boot.subj = TRUE,
                           calc.median.surv = FALSE){

  trt.assign <- match.arg(trt.assign)

  if(calc.median.surv) stop("`Calculation of median survival time has not been implemented yet.")

  if(methods::is(sim, "survparamsim_pre_resampled")){
    if(sim$newdata.orig.missing & calc.obs) {
      warning("Original observed data not provided in `surv_param_sim_pre_resampled()` and KM will not be estimated for the observed data. Speficy `calc.obs = FALSE` to avoid this warning.")
      calc.obs = FALSE
    }
  }

  if(length(trt) > 1) stop("`trt` can only take one string")

  # This needs to be kept as syms - rlang::sym() fails with trt=NULL
  trt.syms   <- rlang::syms(trt)
  group.syms <- rlang::syms(group)
  # This is needed to handle when the same variable is used for both `group` and `trt`
  group.trt.syms <- rlang::syms(unique(c(group, trt)))

  if(length(trt.syms) + length(group.syms) > length(group.trt.syms)){
    warning(paste("Use of the same variable for `group` and `trt` is discouraged.",
                  "If you need a colored & faceted plot, please consider assigning",
                  "your variable to `trt`, and on the plot generated from `plot_km_pi()`,",
                  "apply `facet_wrap()` or `facet_grid()`"))
  }

  ## time for output
  if(is.null(simtimelast)){
    t.out <- seq(0, sim$t.last.orig.new, length.out = 100)
  } else {
    t.out <- seq(0, simtimelast, length.out = round(100 * max(simtimelast/sim$t.last.orig.new, 1)))
  }



  # Observed data --------------------------------

  calc.obs.km.list <- calc_obs_km(sim, calc.obs, group.trt.syms)
  obs.km          <- calc.obs.km.list$obs.km
  obs.median.time <- calc.obs.km.list$obs.median.time


  # Simulated data --------------------------------
  # Calculate percentiles for simulated data

  ## First nest data - KM fit will done for each nested data

  newdata.group <-
    sim$newdata.nona.sim %>%
    dplyr::select(subj.sim, !!!group.trt.syms)

  sim.grouped <-
    sim$sim %>%
    dplyr::left_join(newdata.group, by = "subj.sim") %>%
    dplyr::group_by(rep, !!!group.trt.syms)

  sim.nested <-
    tidyr::nest(sim.grouped) %>%
    dplyr::ungroup()

  # Extract linear predictor (lp), also get scale
  df.lp.extracted <-
    sim.nested %>%
    dplyr::mutate(lp = purrr::map(data, function(x) x$lp)) %>%
    dplyr::select(-data) %>%
    dplyr::left_join(sim$scale.bs.df, by = "rep")


  sim.km <-
    df.lp.extracted %>%
    dplyr::mutate(survfun =
                    purrr::map2(lp, scale,
                                function(x, y) create_survfun(lpvec = x, scale = y, dist = dist))) %>%
    dplyr::select(-lp, -scale) %>%
    dplyr::mutate(km = purrr::map(survfun, function(x) data.frame(time = t.out,
                                                                  surv = x(t.out)))) %>%
    dplyr::arrange(rep, !!!group.trt.syms)



  ## Calc quantile for survival curves
  sim.km.quantile <-
    sim.km %>%
    dplyr::select(rep, !!!group.trt.syms, km) %>%
    tidyr::unnest(km) %>%
    dplyr::group_by(!!!group.trt.syms, time) %>%
    tidyr::nest() %>%
    dplyr::mutate(quantiles = purrr::map(data, function(x)
      dplyr::summarize(x,
                       pi_low = stats::quantile(surv, probs = 0.5 - pi.range/2),
                       pi_med = stats::quantile(surv, probs = 0.5),
                       pi_high= stats::quantile(surv, probs = 0.5 + pi.range/2)))) %>%
    tidyr::unnest(quantiles) %>%
    dplyr::ungroup() %>%
    dplyr::select(-data)



  # Output
  out <- list()

  out$calc.obs <- calc.obs
  out$pi.range   <- pi.range

  out$group.syms <- group.syms
  out$trt.syms    <- trt.syms
  out$group.trt.syms <- group.trt.syms
  out$trt.assign <- trt.assign

  out$simtimelast <- simtimelast
  out$t.last <- sim$t.last.orig.new

  out$obs.km <- obs.km

  out$sim.km <- sim.km
  out$sim.km.quantile <- sim.km.quantile

  structure(out, class = c("survparamsim.kmpi.aveHR", "survparamsim.kmpi"))
}


calc_obs_km <- function(sim, calc.obs, group.trt.syms) {

  out <- list()
  out$obs.km <- NULL
  out$obs.median.time <- NULL

  if(calc.obs){
    # Fit K-M curve to observed data
    obs.grouped <-
      sim$newdata.nona.obs %>%
      dplyr::group_by(!!!group.trt.syms)

    if(length(dplyr::group_vars(obs.grouped)) == 0 &
       utils::packageVersion("tidyr") >= '1.0.0') {
      obs.nested <-
        obs.grouped %>%
        nest2(data = dplyr::everything())
    } else {
      obs.nested <- nest2(obs.grouped)
    }


    ## Define formula
    formula <-
      paste(attributes(formula(sim$survreg))$variables,"~1")[2] %>%
      stats::as.formula()


    ## Calc median and KM curve
    obs.km.nested <-
      obs.nested %>%
      dplyr::mutate(kmfit = purrr::map(data, function(x) survival::survfit(formula, data=x))) %>%
      dplyr::mutate(median = purrr::map_dbl(kmfit, function(x) summary(x)$table["median"]),
                    n      = purrr::map_dbl(kmfit, function(x) summary(x)$table["records"]),
                    km = purrr::map(kmfit, extract_km_obs_time_profile))

    out$obs.km <-
      obs.km.nested %>%
      dplyr::select(-data, -kmfit) %>%
      tidyr::unnest(km) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(surv)) %>%
      dplyr::select(-median)

    out$obs.median.time <-
      obs.km.nested %>%
      dplyr::select(!!!group.trt.syms, median, n)

  }

  return(out)
}




