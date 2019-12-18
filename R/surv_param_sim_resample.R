#' @rdname survparamsim
#' @export
#' @param n.resample Number of subjects for resampled simulations.
#' If `strat.resample` is provided, this needs to be a vector of the length
#' equal to the number of categories in the stratification variable.
#' @param strat.resample String specifying stratification variable for
#' resampling.
surv_param_sim_resample <- function(object, newdata, n.rep = 1000, censor.dur = NULL,
                                    n.resample, strat.resample = NULL){

  # Replace with packageVersion("tidyr") == '1.0.0' if nest issue is resolved in the next version
  # See https://github.com/tidyverse/tidyr/issues/751
  nest2 <- ifelse(utils::packageVersion("tidyr") >= '1.0.0', tidyr::nest_legacy, tidyr::nest)
  unnest2 <- ifelse(utils::packageVersion("tidyr") >= '1.0.0', tidyr::unnest_legacy, tidyr::unnest)

  resample_per_strat <- function(data, n.resample, n.rep){
    dplyr::sample_n(data, n.resample * n.rep, replace = TRUE) %>%
      dplyr::mutate(rep = rep(1:n.rep, each = n.resample))
  }

  if(is.null(strat.resample)){
    newdata.resampled <- resample_per_strat(newdata, n.resample = n.resample, n.rep = n.rep)
  } else {
    strat.sym   <- rlang::sym(strat.resample)

    strat.category <-
      newdata %>%
      dplyr::arrange(!!strat.sym) %>%
      dplyr::pull(!!strat.sym) %>%
      unique()

    if(length(strat.category) != length(n.resample)) {
      stop("Length of `n.resample` needs to be the same as the categories in the statifying variable")
    }

    sample_scheme <- tibble::tibble(!!strat.sym := strat.category,
                                    n.resample = n.resample)

    newdata.resampled <-
      newdata %>%
      dplyr::group_by(!!strat.sym) %>%
      nest2() %>%
      dplyr::left_join(sample_scheme, by = strat.resample) %>%
      dplyr::mutate(sample = purrr::map2(data, n.resample, resample_per_strat, n.rep = n.rep)) %>%
      unnest2(sample)
  }

  newdata.resampled <-
    newdata.resampled %>%
    dplyr::arrange(rep) %>%
    dplyr::mutate(subj.sim.all = dplyr::row_number())

  newdata.resampled.nested <-
    newdata.resampled %>%
    dplyr::group_by(rep) %>%
    nest2()


  simulate_each <- function(data, object, censor.dur){
    sim.each <- surv_param_sim(object, data, n.rep = 1, censor.dur = censor.dur,  na.warning = FALSE)

    sim.each.sim <-
      sim.each$sim %>%
      dplyr::left_join(dplyr::select(sim.each$newdata.nona.sim, subj.sim, subj.sim.all), by = "subj.sim") %>%
      dplyr::mutate(subj.sim = subj.sim.all) %>%
      dplyr::select(-subj.sim.all, -rep)
  }

  sim <-
    newdata.resampled.nested %>%
    dplyr::mutate(sim = purrr::map(data, simulate_each, object = object, censor.dur = censor.dur)) %>%
    unnest2(sim)


  # Generate newdata.nona.obs from non-resample data
  sim.wo.resample <- surv_param_sim(object, newdata, n.rep = 1, censor.dur = censor.dur)


  # Create a list for output
  out <- list()

  out$survreg <- object
  out$newdata <- newdata
  out$newdata.nona.obs <- sim.wo.resample$newdata.nona.obs
  ## Currently NA is not removed in out$newdata.nona.sim, should be fine for now
  out$newdata.nona.sim <- dplyr::rename(newdata.resampled, subj.sim = subj.sim.all)
  out$sim <- sim
  out$censor.dur <- censor.dur

  structure(out, class = c("survparamsim_resample", "survparamsim"))

}
