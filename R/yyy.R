# Set globalVariables to minimize R CMD check notes


if(getRversion() >= "2.15.1"){
  # General
  utils::globalVariables(c('.', ':=', "data", 'subj.sim', 'subj.sim.all', 'median',
                           'pi_low', 'pi_med', 'pi_high', 'obs',
                           'n', 'description', 'quantile', 'event'))

  # calc_hr
  utils::globalVariables(c('coxfit', 'HR'))

  # calc_km
  utils::globalVariables(c('kmfit', 'km', 'surv', 'time', 'cnsr'))

  # extract sim
  utils::globalVariables(c('n.resample', 'covec50', 'cove0', 'covemaxstr', 'covec50str', 'cove0str',
                           'covemaxfct', 'covec50fct', 'cove0fct', 'Covariates'))

  # replace_prm_names
  utils::globalVariables(c('prmname', 'prmname2', 'index', 'prm', 'level'))
}

