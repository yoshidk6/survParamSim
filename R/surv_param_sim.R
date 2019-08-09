#' Simulation of parametric survival model
#'
#' The main function to generate predicted survival using a model object
#' generated with \code{\link[survival]{survreg}} function.
#'
#' @rdname survparamsim
#' @export
#' @param object A `survreg` class object. Currently accept exponential,
#'   lognormal, weibull, loglogistic, and gaussian distributions.
#' @param newdata A requred data frame for simulation that contain covariates in
#'   the survival model. It is required even if this is the same as the one used
#'   for \code{\link[survival]{survreg}} function.
#'
#'   It also has to contain columns for survival information. These can be used
#'   in \code{\link{plot_km_pi}} and \code{\link{plot_hr_pi}} function as
#'   observed data. Survival information can be dummy data, but time need to be
#'   long enough so that simulated KM plot will be long enough for
#'   \code{\link{plot_km_pi}} to draw simulated survival curves.
#'
#'   Subjects with NA for covariates in `survreg` model will be removed from
#'   the simulation and subsequent plotting of observed data.
#' @param n.rep An interger defining numbers of parametric bootstrap runs
#' @param censor.dur A two elements vector specifying duration of events
#'   censoring. Censoring time will be calculated with uniform distribution
#'   between two numbers. No censoring will be applied if NULL is provided.
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
#' \code{\link{surv_param_sim}} returns simulation using the provided subject
#' in `newdata` as it is, while \code{\link{surv_param_sim_resample}} perform
#' simulation based on resampled subjects from the dataset. The latter allows
#' more flexibility in terms of simulating future trials with different number
#' of subjects.
#'
#' Currently we have not tested whether this function work for a `survreg` model
#' with stratification variables.
#'
#'
#' @examples
#' \dontrun{
#'
#' library(survival)
#'
#' fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)
#'
#' object <- fit.lung
#' n.rep  <-  30
#' newdata <- tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog))
#' censor.dur <- c(200, 1100)
#'
#' sim <- surv_param_sim(object, newdata, n.rep, censor.dur)
#'
#' }
#'
#'
#'
surv_param_sim <- function(object, newdata, n.rep = 1000, censor.dur = NULL, na.warning = TRUE){

  if(missing(newdata)) stop("`newdata` needs to be provided even if the same as the one for `survreg()`")

  ## Need to show error if the dataset contains a column named `subj.sim`


  # Prepare model parameter with bootstrap
  thetaHat <- c(object$coef,log(object$scale)) # point estimates of model parameters
  theta    <- mvtnorm::rmvnorm(n.rep, thetaHat, object$var) # parametric bootstrap of model parameters


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
  logmean <- theta[,-ncol(theta)] %*% t(rdata)

  preds <-
    switch(object$dist,
           gaussian =
             rnorm(n    = length(logmean),
                   mean = logmean,
                   sd   = exp(th[,ncol(th)])),
           lognormal =
             rlnorm(n    = length(logmean),
                    mean = logmean,
                    sd   = exp(theta[,ncol(theta)])),
           weibull =
             rweibull(n     = length(logmean),
                      shape = 1/exp(theta[,ncol(theta)]),
                      scale = exp(logmean)),
           loglogistic = exp(
             rlogis(n = length(logmean),
                        location = logmean,
                        scale = exp(theta[,ncol(theta)]))
             ),
           exponential =
             rexp(n = length(logmean),
                  rate = 1/exp(th[]))
    )


  preds <- matrix(preds, nrow=nrow(logmean))
  event.status <- matrix(rep(1,n.rep*n.subj),ncol=ncol(logmean))


  # Censoring time simulation if censor.dur is not NULL
  if(!is.null(censor.dur)){
    censor.time <-
      runif(n.rep*n.subj, censor.dur[1], censor.dur[2]) %>%
      matrix(nrow=nrow(logmean))

    event.status[preds > censor.time] <- 0

    # Replace predicted survival time with censor time if the latter comes earlier
    preds[event.status == 0] <- censor.time[event.status == 0]
  }

  # Export as data frame
  subj.sim.id <-
    tibble::as_tibble(rdata, rownames = "subj.sim") %>%
    .$subj.sim %>%
    as.integer()

  sim <- data.frame(time = as.vector(preds),
                    event= as.vector(event.status),
                    rep  = rep(1:n.rep,  times=n.subj),
                    subj.sim = rep(subj.sim.id,  each =n.rep))

  newdata <-
    newdata %>%
    tibble::rownames_to_column("subj.sim") %>%
    dplyr::mutate(subj.sim = as.integer(subj.sim))

  ## Only keep subjects used for simulation
  ## Subjects with NA for covariates are excluded by model.frame()
  newdata.nona <-
    newdata %>%
    dplyr::filter(subj.sim %in% subj.sim.id)


  # Create a list for output
  out <- list()

  out$survreg <- object
  out$newdata <- newdata
  out$newdata.nona.obs <- newdata.nona
  out$newdata.nona.sim <- newdata.nona
  out$sim <- sim
  out$censor.dur <- censor.dur

  structure(out, class = c("survparamsim"))
}



