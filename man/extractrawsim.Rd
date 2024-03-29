% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_raw_sim.R
\name{extractrawsim}
\alias{extractrawsim}
\alias{extract_sim}
\alias{extract_hr}
\alias{extract_km_obs}
\alias{extract_medsurv}
\alias{extract_medsurv_delta}
\title{Functions to extract raw simulated samples}
\usage{
extract_sim(sim)

extract_hr(hr.pi)

extract_km_obs(km.pi)

extract_medsurv(km.pi)

extract_medsurv_delta(km.pi)
}
\arguments{
\item{sim}{A \code{survparamsim} class object generated by \code{\link[=surv_param_sim]{surv_param_sim()}} function.}

\item{hr.pi}{a return object from \code{\link[=calc_hr_pi]{calc_hr_pi()}} function.}

\item{km.pi}{A return object from \code{\link[=calc_km_pi]{calc_km_pi()}} function.}
}
\description{
Functions to extract raw simulated samples
}
\details{
\code{\link[=extract_sim]{extract_sim()}} extracts raw survival time & event status for all simulated subjects.

\code{\link[=extract_hr]{extract_hr()}} extracts simulated HR for all repeated simulations.
If HR was calculated based on simulated survival times with
\code{\link[=calc_hr_pi]{calc_hr_pi()}} function, it also returns p values for Cox regression
fits, one for each group
based on Wald test and another for the overall significance of the
coefficient based on logrank test. The latter has the same values
across treatment groups when >2 levels in treatment

\code{\link[=extract_km_obs]{extract_km_obs()}} extracts observed Kaplan-Meier curves.

\code{\link[=extract_medsurv]{extract_medsurv()}} extracts simulated median survival times for all repeated simulations

\code{\link[=extract_medsurv_delta]{extract_medsurv_delta()}} extracts delta of median survival times between treatment groups
}
