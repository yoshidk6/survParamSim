

# survParamSim 0.1.5 (in development)

## Breaking changes

The original `extract_median_surv` function was renamed to `extract_median_surv_pi` for consistent name conventions  
The new `extract_median_surv` function will return raw median survival time per simulation instead of prediction interval  

## Major changes

The new `extract_median_surv_delta_pi` function will return prediction intervals of the difference of median survival time between treatment  

## Minor changes


# survParamSim 0.1.4

## Minor changes

* Bug fix - use `vdiffr` conditionally on tests


# survParamSim 0.1.3 

## Minor changes

* Bug fix - disable `options(lifecycle_verbosity = "error")` to avoid unnecessary errors


# survParamSim 0.1.2

## Minor changes

* Bug fix of `surv_param_sim_resample()` with `tidyr >= 1.0.0`

# survParamSim 0.1.1

## Minor changes

* Fix unit test code to work with dplyr 1.0.0

* Add example outputs in Readme and vignette


# survParamSim 0.1.0

* Initial release
