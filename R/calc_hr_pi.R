#' Generate hazard ratio with prediction intervals from parametric bootstrap simulation
#'
#' @name calculate_hr_pi
NULL

#' @rdname calculate_hr_pi
#' @export
#' @param sim A `survparamsim` class object generated by [surv_param_sim()] function.
#' @param trt A string to specify which column define treatment status to calculate HR.
#' @param group Optional string(s) to specify grouping variable(s).
#' You will have faceted histograms for these variables in [plot_hr_pi()] function.
#' @param pi.range Prediction interval for simulated HR.
#' @param calc.obs A logical to specify whether to calculate HR for the observed data.
#' Need be set as FALSE if survival information in the `newdata` is dummy.
#' @param trt.assign Specify which of the categories of `trt` need to be considered as control group.
#' See details below if you have more than two categories.
#'
#' @details
#' [calc_hr_pi()] calculate hazard ratio using the simulated survival times with Cox proportional hazard
#' model, while [calc_ave_hr_pi()] calculate "average" hazard ratio using the mean survival & probability
#' density function per treatment groups.
#'
#' If your `trt` has more than two categories/levels and want to specify which one to use as a
#' reference group, you can convert the column into a factor in the `newdata` input for
#' [surv_param_sim()]. The first level will be used as a reference group.
#'
calc_hr_pi <- function(sim, trt, group = NULL, pi.range = 0.95,
                       calc.obs = TRUE, trt.assign = c("default", "reverse")){

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
                            calc.obs)


  # Calc HR for simulated data ----------------------------------------------------------------

  ## First nest data - cox fit will done for each nested data
  newdata.trt.group <-
    newdata.nona.sim %>%
    dplyr::select(subj.sim, !!trt.sym, !!!group.syms)

  sim.nested <-
    sim$sim %>%
    dplyr::left_join(newdata.trt.group, by = "subj.sim") %>%
    dplyr::group_by(rep, !!!group.syms) %>%
    tidyr::nest()


  # Calculate HR with Cox
  sim.hr <- calc_hr_for_sim_with_cox(sim.nested, trt, trt.sym, trt.levels)


  ## Reverse back the factor
  if(trt.assign == "reverse"){
    sim.hr <-
      sim.hr %>%
      dplyr::mutate(!!trt.sym := forcats::fct_rev(!!trt.sym))
  }


  # Calc quantiles ----------------------------------------------------------------

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

  structure(out, class = c("survparamsim.hrpi.simulated.time", "survparamsim.hrpi"))
}


#' Plot simulated HR histogram(s) overlayed with prediction intervals
#'
#' @export
#' @param hr.pi a return object from \code{\link{calc_hr_pi}} function.
#' @param show.obs A logical specifying whether to show observed HR on the plot.
#'   This will have no effect if `calc.obs` was set to `FALSE` in \code{\link{calc_hr_pi}}.
#'
plot_hr_pi <- function(hr.pi, show.obs = TRUE){

  obs.hr <- hr.pi$obs.hr
  sim.hr <- hr.pi$sim.hr
  hr.pi.quantile  <- hr.pi$hr.pi.quantile
  trt.levels <- hr.pi$trt.levels

  group.syms <- hr.pi$group.syms
  trt.sym    <- hr.pi$trt.sym

  g <-
    ggplot2::ggplot(sim.hr, ggplot2::aes(HR)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.5, color = "black") +
    ggplot2::geom_vline(data = dplyr::filter(hr.pi.quantile, description %in% c("pi_low", "pi_high")),
                        ggplot2::aes(xintercept = HR),
                        lty="dashed")


  ## Observed
  if(hr.pi$calc.obs & show.obs) {
    g <-
      g +
      ggplot2::geom_vline(data = dplyr::filter(hr.pi.quantile, description == "obs"),
                          ggplot2::aes(xintercept = HR),
                          color = "red", lwd = 1)

  }

  # Facet fig based on group
  if(length(trt.levels) > 2) {
    g <- g + ggplot2::facet_grid(ggplot2::vars(!!trt.sym),
                                 ggplot2::vars(!!!group.syms),
                                 labeller = ggplot2::label_both)
  } else {
    if(length(group.syms) == 1 || length(group.syms) >= 3 ) {
      g <- g + ggplot2::facet_wrap(ggplot2::vars(!!!group.syms),
                                   labeller = ggplot2::label_both)
    } else if (length(group.syms) == 2) {
      g <- g + ggplot2::facet_grid(ggplot2::vars(!!group.syms[[1]]),
                                   ggplot2::vars(!!group.syms[[2]]),
                                   labeller = ggplot2::label_both)
    }
  }

  return(g)

}


handle_trt_group <- function(trt, trt.assign, group, newdata.nona.obs, newdata.nona.sim){
  if(missing(trt)) stop("`trt` needs to be specified")
  if(length(trt) > 1) stop("`trt` can only take one string")

  group.syms <- rlang::syms(group)
  trt.sym    <- rlang::sym(trt)

  # Check trt values
  n.distinct.trt <- check_trt(newdata.nona.obs, trt.sym)

  # Convert trt to factor
  newdata.nona.obs <-
    newdata.nona.obs %>%
    dplyr::mutate(!!trt.sym := factor(!!trt.sym))
  newdata.nona.sim <-
    newdata.nona.sim %>%
    dplyr::mutate(!!trt.sym := factor(!!trt.sym))

  # Reverse control vs trt
  if(trt.assign == "reverse"){
    newdata.nona.obs <-
      newdata.nona.obs %>%
      dplyr::mutate(!!trt.sym := forcats::fct_rev(!!trt.sym))
    newdata.nona.sim <-
      newdata.nona.sim %>%
      dplyr::mutate(!!trt.sym := forcats::fct_rev(!!trt.sym))
  }

  trt.levels <- dplyr::pull(newdata.nona.obs, !!trt.sym) %>% levels()

  out <- list()
  out$group.syms       <- group.syms
  out$trt.sym          <- trt.sym
  out$newdata.nona.obs <- newdata.nona.obs
  out$newdata.nona.sim <- newdata.nona.sim
  out$trt.levels       <- trt.levels

  return(out)
}


check_trt <- function(newdata.nona.obs, trt.sym){

  # Check trt values
  trt.vec <-
    newdata.nona.obs %>%
    dplyr::select(!!trt.sym) %>%
    .[[1]]

  n.distinct.trt <- length(unique(trt.vec))

  if(sum(is.na(trt.vec)) > 0) stop("`trt` cannot has NA values")
  if(is.factor(trt.vec) & nlevels(trt.vec) != n.distinct.trt) warning("`trt` variable is factor and has unused levels, which is automatically dropped`")
  if(is.ordered(trt.vec)) stop("`trt` cannot be an ordered factor, please use a regular factor instead")

  return(n.distinct.trt)

}


calc_hr_for_obs <- function(sim, newdata.nona.obs, group.syms, trt, trt.sym, trt.assign, trt.levels,
                            calc.obs) {

  if(calc.obs){

    obs.grouped <-
      newdata.nona.obs %>%
      dplyr::group_by(!!!group.syms)

    if(length(dplyr::group_vars(obs.grouped)) == 0) {
      obs.nested <-
        obs.grouped %>%
        tidyr::nest(data = dplyr::everything())
    } else {
      obs.nested <- tidyr::nest(obs.grouped)
    }

    ## Define function to calc HR
    calc_hr_each_obs <- function(x){
      formula <-
        paste(attributes(formula(sim$survreg))$variables,"~",trt)[2] %>%
        stats::as.formula()

      cfit <- survival::coxph(formula, data=x)
      p.value.logrank <- broom::glance(cfit)$p.value.log
      cfit %>%
        broom::tidy(exponentiate = TRUE) %>%
        dplyr::mutate(!!trt.sym := factor(substr(term, nchar(trt)+1, nchar(term)),
                                          levels = trt.levels)) %>%
        dplyr::select(!!trt.sym, HR = estimate, p.value.coef.wald = p.value) %>%
        dplyr::mutate(p.value.logrank = p.value.logrank)
    }
    safe_calc_hr_each_obs <- purrr::safely(calc_hr_each_obs, otherwise = NA)

    ## Calc HR
    obs.hr <-
      obs.nested %>%
      dplyr::mutate(coxfit = purrr::map(data, safe_calc_hr_each_obs),
                    HR = purrr::map(coxfit, ~.$result),
                    description = "obs") %>%
      tidyr::unnest(HR) %>%
      dplyr::select(-data, -coxfit) %>%
      dplyr::ungroup()

    # Reverse back the factor
    if(trt.assign == "reverse"){
      obs.hr <-
        obs.hr %>%
        dplyr::mutate(!!trt.sym := forcats::fct_rev(!!trt.sym))
    }

    if(dplyr::filter(obs.hr, is.na(HR)) %>% nrow() > 0) {
      warning("HR was not calculable in at least one subgroup for the observed data, likely due to small number of subjects.")
    }

  } else {
    obs.hr <- NULL
  }

  return(obs.hr)
}


calc_hr_for_sim_with_cox <- function(sim.nested, trt, trt.sym, trt.levels) {

  ## Define function to calc HR
  calc_hr_each_sim <- function(x){
    formula <-
      paste("Surv(time, event) ~",trt) %>%
      stats::as.formula()

    cfit <- survival::coxph(formula, data=x)
    p.value.logrank <- broom::glance(cfit)$p.value.log
    cfit %>%
      broom::tidy(exponentiate = TRUE) %>%
      dplyr::mutate(!!trt.sym := factor(substr(term, nchar(trt)+1, nchar(term)),
                                        levels = trt.levels))%>%
      dplyr::select(!!trt.sym, HR = estimate, p.value.coef.wald = p.value) %>%
      dplyr::mutate(p.value.logrank = p.value.logrank)
  }
  safe_calc_hr_each_sim <- purrr::safely(calc_hr_each_sim, otherwise = NA)

  ## Calc HR
  sim.hr <-
    sim.nested %>%
    dplyr::mutate(coxfit = purrr::map(data, safe_calc_hr_each_sim),
                  HR = purrr::map(coxfit, ~.$result),
                  description = "sim") %>%
    tidyr::unnest(HR) %>%
    dplyr::select(-data, -coxfit) %>%
    dplyr::ungroup()

  if(dplyr::filter(sim.hr, is.na(HR)) %>% nrow() > 0) {
    warning("HR was not calculable in at least one subgroup for the simulated data, likely due to small number of subjects and this might results in biased estimates for prediction intervals.")
  }

  return(sim.hr)
}


calc_hr_quantiles <- function(pi.range, sim.hr, obs.hr, calc.obs,
                              group.syms, trt.sym) {

  quantiles <-
    tibble::tibble(description = c("pi_low", "pi_med", "pi_high"),
                   quantile = c(0.5 - pi.range/2, 0.5, 0.5 + pi.range/2))

  sim.hr.pi <-
    sim.hr %>%
    dplyr::group_by(!!!group.syms, !!trt.sym) %>%
    dplyr::summarize(pi_low = as.numeric(stats::quantile(HR, probs = 0.5 - pi.range/2, na.rm = TRUE)),
                     pi_med = as.numeric(stats::quantile(HR, probs = 0.5, na.rm = TRUE)),
                     pi_high= as.numeric(stats::quantile(HR, probs = 0.5 + pi.range/2, na.rm = TRUE))) %>%
    dplyr::ungroup() %>%
    tidyr::gather(description, HR, pi_low:pi_high) %>%
    dplyr::left_join(quantiles, by = "description")

  if(calc.obs){
    hr.pi.quantile <-
      dplyr::select(obs.hr, -p.value.coef.wald, -p.value.logrank) %>%
      dplyr::bind_rows(sim.hr.pi, .) %>%
      dplyr::arrange(!!!group.syms, !!trt.sym)

  } else {
    hr.pi.quantile <-
      sim.hr.pi %>%
      dplyr::arrange(!!!group.syms, !!trt.sym)
  }

  return(hr.pi.quantile)
}


#' @rdname survparamsim-methods
#' @export
print.survparamsim.hrpi <- function(x, ...){
  trt <- as.character(x$trt.sym)
  group <- as.character(x$group.syms)

  cat("---- Simulated and observed (if calculated) hazard ratio ----\n")
  cat("* Use `extract_hr_pi()` to extract prediction intervals and observed HR\n")
  cat("* Use `extract_hr()` to extract individual simulated HRs\n")
  cat("* Use `plot_hr_pi()` to draw histogram of predicted HR\n\n")
  cat("* Settings:\n")
  cat("    trt: ", trt, "\n", sep="")
  cat("         ('", paste(x$trt.levels[-1], collapse = "', '"), "' as test trt, '", x$trt.levels[[1]], "' as control)\n", sep="")
  cat("    group:", ifelse(is.null(group), "(NULL)", group), "\n", sep=" ")
  cat("    pi.range:", x$pi.range, "\n", sep=" ")
  cat("    calc.obs:", x$calc.obs, "\n", sep=" ")

  if(methods::is(x, "survparamsim.hrpi.simulated.time")){
    cat("    method: Cox on simulated survival using `calc_hr_pi()`")
  } else if (methods::is(x, "survparamsim.hrpi.aveHR")){
    cat("    method: average HR using `calc_ave_hr_pi()`")
  }

}


#' @rdname survparamsim-methods
#' @export
summary.survparamsim.hrpi <- function(object, ...) {

  return(extract_hr_pi(object))
}



