
# survParamSim 0.1.6

## Minor changes

* Fix error when the same variable is included in both `group` and `trt` for `calc_km_pi()`

# survParamSim 0.1.5

## Minor changes

* New function `extract_medsurv_delta_pi()` will return prediction intervals of the difference of median survival time between treatment  

* The `extract_median_surv()` function was deprecated and superseded by `extract_medsurv_pi()` for consistent name conventions  

* New function `extract_medsurv()` will return raw median survival time per simulation

* New function `surv_param_sim_pre_resampled()` enables simulation with an already resampled dataset

* New function `calc_km_pi()` shows warning when there are survival simulations for which median survival times were not reached

* `calc_hr_pi()` can now calculate hazard ratios for dataset with more than two treatment groups


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
