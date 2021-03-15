#' @rdname survparamsim_pre_resampled
#' @export
#' @param object A `survreg` class object. Currently accept exponential,
#'   lognormal, weibull, loglogistic, and gaussian distributions.
#' @param newdata.orig AAAA
#'
#'   Subjects with NA for covariates in `survreg` model are not allowed for this
#'   function.
#' @param newdata.resampled AAAA
#'   already resampled
#'   Some specifications:
#'     - rep must be used
#'     - Same number of subjects per rep (if subgroup in HR or KM, stratify for these variables too)
#' @param censor.dur A two elements vector specifying duration of events
#'   censoring. Censoring time will be calculated with uniform distribution
#'   between two numbers. No censoring will be applied if NULL is provided.
#' @return A `survparamsim` object that contains the original `survreg` class
#'   object, newdata, and a data frame for predicted survival profiles.
#' @details
#' See \code{\link{surv_param_sim}} for additional details.
#' \code{\link{surv_param_sim_pre_resampled}} performs simulations on an already-
#'  resampled dataset.
#'
#'
surv_param_sim_pre_resampled <- function(object, newdata.orig, newdata.resampled, censor.dur = NULL){




  ###############
  ###############
  # test_that("Check N consistent across reps", {
  # })
  #
  check_data_n_per_resample(newdata.resampled)

  check_data_na_resample(newdata.orig, object)
  check_data_na_resample(newdata.resampled, object)

  ###############
  ###############

  newdata.resampled <-
    newdata.resampled %>%
    dplyr::arrange(rep) %>%
    dplyr::rename(.rep.pre.resample = rep) %>%
    dplyr::mutate(subj.sim.all = dplyr::row_number())

  sim <-
    surv_param_sim(object, newdata = newdata.resampled, n.rep = 1, censor.dur = censor.dur,  na.warning = TRUE)

  sim.sim <-
    sim$sim %>%
    dplyr::left_join(dplyr::select(extract_sim(sim), subj.sim, .rep.pre.resample, subj.sim.all), by = "subj.sim") %>%
    dplyr::select(-rep, -subj.sim) %>%
    dplyr::rename(rep = .rep.pre.resample,
                  subj.sim = subj.sim.all)


  # Generate newdata.nona.obs and t.last.orig.new from non-resample data
  sim.wo.resample <- surv_param_sim(object, newdata.orig, n.rep = 1, censor.dur = censor.dur)


  # Create a list for output
  out <- list()

  out$survreg <- object
  out$t.last.orig.new <- sim.wo.resample$t.last.orig.new
  out$newdata.nona.obs <- sim.wo.resample$newdata.nona.obs
  out$newdata.nona.sim <- dplyr::rename(newdata.resampled,
                                        rep = .rep.pre.resample,
                                        subj.sim = subj.sim.all) # Used for grouping assignment in sim HR & K-M
  out$sim <- sim.sim
  out$censor.dur <- censor.dur

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

