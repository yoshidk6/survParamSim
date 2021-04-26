#' Simulation of parametric survival model
#'
#' The main function(s) to generate predicted survival using a model object
#' generated with [survival::survreg()] function.
#'
#' @name survparamsim
NULL


#' @rdname survparamsim
#' @export
#' @param object A `survreg` class object. Currently accept exponential,
#'   lognormal, weibull, loglogistic, and gaussian distributions.
#' @param newdata A required data frame for simulation that contain covariates in
#'   the survival model. It is required even if this is the same as the one used
#'   for [survival::survreg] function.
#'
#'   It also has to contain columns for survival information. These can be used
#'   in [plot_km_pi()] and [plot_hr_pi()] function as
#'   observed data. Survival information can be dummy data, but time need to be
#'   long enough so that simulated KM plot will be long enough for
#'   [plot_km_pi()] to draw simulated survival curves.
#'
#'   Subjects with NA for covariates in `survreg` model will be removed from
#'   the simulation and subsequent plotting of observed data.
#' @param n.rep An integer defining numbers of parametric bootstrap runs
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
#'   object, newdata, and a data frame for predicted survival profiles with the
#'   following columns:
#'   \itemize{\item \strong{time}: predicted event or censor time
#'            \item \strong{event}: event status, 0=censored, 1=event
#'            \item \strong{rep}: ID for parametric bootstrap runs
#'            \item \strong{subj}: ID for subjects in newdata (currently
#'              original ID is not retained and subj is sequentially assigned
#'              as 1:nrow(newdata))
#'   }
#'
#' @details
#' [surv_param_sim()] returns simulation using the provided subject
#' in `newdata` as it is, while \code{\link{surv_param_sim_resample}} perform
#' simulation based on resampled subjects from the dataset. The latter allows
#' more flexibility in terms of simulating future trials with different number
#' of subjects.
#' Note that with [surv_param_sim_resample()], there is no
#' automatic safeguard to ensure certain number of subjects in each subgroup
#' or treatment groups, which may result in inconsistent number of subjects per
#' simulation or leads to Cox regression instability due to small N. Consider
#' using stratified resampling in this case.
#'
#' Currently we have not tested whether this function work for a `survreg` model
#' with stratification variables.
#'
#'
#' @examples
#' library(survival)
#'
#' fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)
#'
#' object <- fit.lung
#' n.rep  <-  30
#' newdata <-
#'   tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
#'   tidyr::drop_na()
#' censor.dur <- c(200, 1100)
#'
#' sim <- surv_param_sim(object, newdata, n.rep, censor.dur)
#'
#'
#'
surv_param_sim <- function(object, newdata, n.rep = 1000, censor.dur = NULL,
                           coef.var = TRUE, na.warning = TRUE){

  if(missing(newdata)) stop("`newdata` needs to be provided even if the same as the one for `survreg()`")

  ## Need to show error if the dataset contains a column named `subj.sim`


  check_censor_dur(censor.dur)

  # Prepare model parameter with bootstrap
  ## point estimates of model parameters
  if(object$dist != "exponential"){
    theta <- c(object$coef, log(object$scale))
  } else {
    theta <- object$coef
  }

  # parametric bootstrap of model parameters
  if(coef.var) {
    th.bs <- mvtnorm::rmvnorm(n.rep, theta, stats::vcov(object))
  } else {
    th.bs <- matrix(rep(theta, each = n.rep), nrow = n.rep)
  }


  # Prep data matrix for simulation
  mf <- stats::model.frame(object, data = newdata)
  rdata <- stats::model.matrix(object, mf)
  ## Get #subjects and subject IDs, after NA excluded by model.frame()
  n.subj <- nrow(mf)

  if(n.subj == 0) {
    stop("No subjects present in `newdata` for simulation. It might be because all subjects has NA in model variables (including survival and censoring status)")
  }
  if(n.subj < nrow(newdata) & na.warning) {
    warning("Not all subjects in `newdata` will be used for `surv_param_sim`, likely because NA is present")
  }


  # Survival time simulation
  ## Generate linear predictors (lp)
  if(object$dist != "exponential"){
    lp <- th.bs[,-ncol(th.bs)] %*% t(rdata)
    scale.bs <- th.bs[,ncol(th.bs)]
  } else {
    lp <- th.bs %*% t(rdata)
  }

  preds <-
    switch(object$dist,
           gaussian = stats::rnorm(n    = length(lp),
                                   mean = lp,
                                   sd   = exp(scale.bs)),
           lognormal = stats::rlnorm(n    = length(lp),
                                     mean = lp,
                                     sd   = exp(scale.bs)),
           weibull = stats::rweibull(n     = length(lp),
                                     shape = 1/exp(scale.bs),
                                     scale = exp(lp)),
           loglogistic = exp(stats::rlogis(n = length(lp),
                                           location = lp,
                                           scale = exp(scale.bs))),
           exponential = stats::rexp(n = length(lp),
                                     rate = 1/exp(lp))
    )


  preds <- matrix(preds, nrow=nrow(lp))
  event.status <- matrix(rep(1,n.rep*n.subj),ncol=ncol(lp))


  # Censoring time simulation if censor.dur is not NULL
  if(!is.null(censor.dur)){
    censor.time <-
      stats::runif(n.rep*n.subj, censor.dur[1], censor.dur[2]) %>%
      matrix(nrow = nrow(lp))

    event.status[preds > censor.time] <- 0

    # Replace predicted survival time with censor time if the latter comes earlier
    preds[event.status == 0] <- censor.time[event.status == 0]
  }

  # Export as data frame
  subj.sim.id <-
    tibble::as_tibble(rdata, rownames = "subj.sim") %>%
    .$subj.sim %>%
    as.integer()

  sim <-
    data.frame(time = as.vector(preds),
               event= as.vector(event.status),
               rep  = rep(1:n.rep,  times=n.subj),
               subj.sim = rep(subj.sim.id,  each =n.rep)) %>%
    dplyr::arrange(rep, subj.sim)

  newdata <-
    newdata %>%
    tibble::rownames_to_column("subj.sim") %>%
    dplyr::mutate(subj.sim = as.integer(subj.sim))

  ## Only keep subjects used for simulation
  ## Subjects with NA for covariates are excluded by model.frame()
  newdata.nona <-
    newdata %>%
    dplyr::filter(subj.sim %in% subj.sim.id)


  ## Last observed time to be used for KM calculations
  ### Calculate last time from original dataset just in case newdata's survival data is dummy
  formula <-
    paste(attributes(formula(object))$variables,"~1")[2] %>%
    stats::as.formula()

  t.last.newdata  <- survival::survfit(formula, data = newdata) %>% .$time
  t.last.origdata <- as.numeric(object$y[,1])

  t.last.orig.new <-max(c(t.last.newdata, t.last.origdata))

  # Create a list for output
  out <- list()

  out$survreg <- object
  out$t.last.orig.new <- t.last.orig.new
  out$newdata.nona.obs <- newdata.nona # Used for obs HR & KM calculation and check TRT variables in calc_hr_pi
  out$newdata.nona.sim <- newdata.nona # Used for grouping assignment in sim HR & K-M
  out$sim <- sim
  out$censor.dur <- censor.dur

  structure(out, class = c("survparamsim"))
}


check_censor_dur <- function(censor.dur = NULL) {
  if(!is.null(censor.dur)) {
    if(length(censor.dur) != 2) stop("censor.dur has to be length two vector or NULL")
    if(censor.dur[1] > censor.dur[2]) stop("censor.dur[2] has to be larger than censor.dur[1]")
  }
}

