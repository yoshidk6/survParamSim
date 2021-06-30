context("test-surv_param_sim")

library(survival)
set.seed(12345)

fit.lung <- survreg(Surv(time, status) ~ sex + ph.ecog, data = lung)

object <- fit.lung
n.rep  <-  30
newdata <-
  tibble::as_tibble(dplyr::select(lung, time, status, sex, ph.ecog)) %>%
  dplyr::mutate(.id.orig = dplyr::row_number())
censor.dur <- c(200, 1100)


sim <- suppressWarnings(surv_param_sim(object, newdata, n.rep, censor.dur))

# sim <- surv_param_sim(object, newdata, n.rep = 1, censor.dur)


test_that("have NA in dataset", {
  expect_warning(surv_param_sim(object, newdata, n.rep, censor.dur),
                 "Not all subjects in `newdata`")
})



test_that("make sure missing NA subjects were removed", {
  ids.in.sim <- extract_sim(sim) %>% dplyr::pull(.id.orig) %>% unique()
  expect_false(14 %in% ids.in.sim) # 14 is subject with NA in ph.ecog
})


# Check multiple imputation samples can run w/o error
set.seed(12345)

## Prep data
m.imp <- 5

data.mi <-
  newdata %>%
  dplyr::select(time, status, sex, ph.ecog) %>%
  dplyr::mutate(time = log(time))

predictormat <- 1 - diag(1,ncol(data.mi));
predictormat[,2] <- 0;

imp.1 <- mice::mice(dplyr::filter(data.mi, status == 1), m = m.imp, predictorMatrix=predictormat, seed=100, printFlag = FALSE);
imp.2 <- mice::mice(dplyr::filter(data.mi, status == 2), m = m.imp, predictorMatrix=predictormat, seed=100, printFlag = FALSE);

imp <- mice::rbind(imp.1, imp.2)

## Fit and put together summary
mit <- mitools::imputationList(lapply(1:m.imp, mice::complete, data=imp))
models <- with(mit, survreg(Surv(exp(time), status) ~ sex + ph.ecog))
betas <- mitools::MIextract(models,fun=coef)
for(i in 1:m.imp) betas[[i]] <- c(betas[[i]], `Log(scale)`=log(models[[i]]$scale)) # Add log(scale) parameter at the end. This needs to be ignored for exponential model
vars <- mitools::MIextract(models, fun=vcov)
fit.mi <- mitools::MIcombine(betas,vars)

sim <- suppressWarnings(surv_param_sim(models[[1]], newdata, n.rep, censor.dur, mi.resuls = fit.mi))


