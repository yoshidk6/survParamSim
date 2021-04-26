#' Simulation of parametric survival model with an already-resampled dataset
#'
#' @export
#' @param object A `survreg` class object. Currently accept exponential,
#'   lognormal, weibull, loglogistic, and gaussian distributions.
#' @param newdata.resampled A required input, the already resampled dataset for simulation.
#'  This dataset must have: (a) `rep` variable indicating the #simulation groups, and (b) the same number of subjects per each `rep`
#' @param newdata.orig  An optional input needed for calculating KM and HR for the observed data.
#' @param censor.dur A two elements vector specifying duration of events
#'   censoring. Censoring time will be calculated with uniform distribution
#'   between two numbers. No censoring will be applied if NULL is provided.
#' @param coef.var Boolean specifying whether parametric bootstrap are
#' performed on survival model coefficients, based on variance-covariance
#' matrix. If FALSE, prediction interval only reflects inherent variability
#' from survival events.
#' @param na.warning Boolean specifying whether warning will be shown if
#' `newdata` contain subjects with missing model variables.
#' @return A `survparamsim` object that contains the original `survreg` class
#'   object, newdata, and a data frame for predicted survival profiles.
#' @details
#' See [surv_param_sim()] for additional details.
#'
#'
surv_param_sim_pre_resampled <- function(object, newdata.resampled, newdata.orig = NULL, censor.dur = NULL,
                                         coef.var = TRUE, na.warning = TRUE){

  # Replace nest with packageVersion("tidyr") == '1.0.0' for a speed issue
  # See https://github.com/tidyverse/tidyr/issues/751
  nest2 <- ifelse(utils::packageVersion("tidyr") == '1.0.0', tidyr::nest_legacy, tidyr::nest)
  unnest2 <- ifelse(utils::packageVersion("tidyr") == '1.0.0', tidyr::unnest_legacy, tidyr::unnest)

  if(is.null(newdata.orig)){
    newdata.orig <- newdata.resampled
    newdata.orig.missing <- TRUE
  } else {
    newdata.orig.missing <- FALSE
  }

  check_data_n_per_resample(newdata.resampled)

  # check_data_na_resample(newdata.orig, object)
  # check_data_na_resample(newdata.resampled, object)



  newdata.resampled <-
    newdata.resampled %>%
    dplyr::arrange(rep) %>%
    dplyr::mutate(subj.sim.all = dplyr::row_number())


  newdata.resampled.nested <-
    newdata.resampled %>%
    dplyr::group_by(rep) %>%
    nest2()


  simulate_each <- function(data, object, censor.dur){
    sim.each <- surv_param_sim(object, data, n.rep = 1, censor.dur = censor.dur,
                               coef.var = coef.var, na.warning = FALSE)

    sim.each.sim <-
      sim.each$sim %>%
      dplyr::left_join(dplyr::select(sim.each$newdata.nona.sim, subj.sim, subj.sim.all), by = "subj.sim") %>%
      dplyr::mutate(subj.sim = subj.sim.all) %>%
      dplyr::select(-subj.sim.all, -rep)
  }

  sim <-
    newdata.resampled.nested %>%
    dplyr::mutate(sim = purrr::map(data, simulate_each, object = object, censor.dur = censor.dur)) %>%
    dplyr::select(-data) %>%
    unnest2(sim) %>%
    dplyr::ungroup()


  # Generate newdata.nona.obs and t.last.orig.new from non-resample data
  sim.wo.resample <- surv_param_sim(object, newdata.orig, n.rep = 1, censor.dur = censor.dur,
                                    coef.var = coef.var, na.warning = na.warning)


  # Create a list for output
  out <- list()

  out$survreg <- object
  out$t.last.orig.new <- sim.wo.resample$t.last.orig.new
  out$newdata.nona.obs <- sim.wo.resample$newdata.nona.obs
  out$newdata.nona.sim <- dplyr::rename(newdata.resampled,
                                        subj.sim = subj.sim.all) # Used for grouping assignment in sim HR & K-M
  out$sim <- sim
  out$censor.dur <- censor.dur

  out$newdata.orig.missing <- newdata.orig.missing

  structure(out, class = c("survparamsim_pre_resampled", "survparamsim"))
}





check_data_n_per_resample <- function(data) {

  n.of.unique.n.per.group <-
    data %>%
    dplyr::count(rep) %>%
    dplyr::pull(n) %>%
    unique() %>%
    length()

  if(n.of.unique.n.per.group > 1) {
    stop("#Subjects must be the same across resampled groups defined by `rep`")
  }

}

