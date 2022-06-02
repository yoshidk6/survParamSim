# Set globalVariables to minimize R CMD check notes


if(getRversion() >= "2.15.1"){
  # General
  utils::globalVariables(c('.', ':=', "data", 'subj.sim', 'subj.sim.all', 'median',
                           'pi_low', 'pi_med', 'pi_high', 'obs',
                           'n', 'description', 'quantile', 'event'))

  # calc_hr
  utils::globalVariables(c('coxfit', 'HR', 'term', 'estimate',
                           "p.value", "p.value.coef.wald", "p.value.logrank"))

  # calc_km
  utils::globalVariables(c('kmfit', 'km', 'surv', 'time', 'cnsr',
                           'is.median.na', 'N.median.NA', 'N.all', 'n_min', 'n_max',
                           '.trt.control.group', '.trt.group.1', '.trt.group.index', 'median_delta'))

  # extract sim
  utils::globalVariables(c('n.resample', 'covec50', 'cove0', 'covemaxstr', 'covec50str', 'cove0str',
                           'covemaxfct', 'covec50fct', 'cove0fct', 'Covariates',
                           '.delta'))

  # replace_prm_names
  utils::globalVariables(c('prmname', 'prmname2', 'index', 'prm', 'level'))
}

