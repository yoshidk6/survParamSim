
# Create a function to calculate S(t)
# The function to be created will take time (`x`) as the input argument and return
create_survfun <- function(lpvec, scale.ln,
                           dist = c("lognormal",
                                    "gaussian",
                                    "weibull",
                                    "loglogistic",
                                    "exponential")){

  dist <- match.arg(dist)

  if(dist == "lognormal") {
    function(x) {
      survmatrix <-
        vapply(lpvec,
               function(lp) {stats::plnorm(q=x, meanlog=lp, sdlog=exp(scale.ln), lower=FALSE)},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else if(dist == "gaussian") {
    function(x) {
      survmatrix <-
        vapply(lpvec,
               function(lp) {stats::pnorm(q=x, mean=lp, sd=exp(scale.ln), lower=FALSE)},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else if(dist == "weibull") {
    function(x) {
      survmatrix <-
        vapply(lpvec,
               function(lp) {stats::pweibull(q=x, shape=1/exp(scale.ln), scale=exp(lp), lower=FALSE)},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else if(dist == "exponential") {
    function(x) {
      survmatrix <-
        vapply(lpvec,
               function(lp) {stats::pexp(q=x, rate=1/exp(lp), lower=FALSE)},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else if(dist == "loglogistic") {
    function(x) {
      survmatrix <-
        vapply(lpvec,
               function(lp) {eha::pllogis(q=x, shape = 1/exp(scale.ln), scale = exp(lp), lower=FALSE)},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else {
    stop("Distribution not defined")
  }
}



# Create a function to calculate PDF(t)
create_pdf <- function(lpvec, scale.ln, dist = "lognormal"){


  if(dist == "lognormal") {
    function(x) {
      pdmatrix <-
        vapply(lpvec,
               function(lp) {stats::dlnorm(x=x, meanlog=lp, sdlog=exp(scale.ln))},
               numeric(length(x)))

      if(length(x) == 1) return(mean(pdmatrix))
      return(rowMeans(pdmatrix))
    }
  } else if(dist == "gaussian") {
    function(x) {
      pdmatrix <-
        vapply(lpvec,
               function(lp) {stats::pnorm(x=x, mean=lp, sd=exp(scale.ln))},
               numeric(length(x)))

      if(length(x) == 1) return(mean(pdmatrix))
      return(rowMeans(pdmatrix))
    }
  } else if(dist == "weibull") {
    function(x) {
      pdmatrix <-
        vapply(lpvec,
               function(lp) {stats::dweibull(x=x, shape=1/exp(scale.ln), scale=exp(lp))},
               numeric(length(x)))

      if(length(x) == 1) return(mean(pdmatrix))
      return(rowMeans(pdmatrix))
    }
  } else if(dist == "exponential") {
    function(x) {
      pdmatrix <-
        vapply(lpvec,
               function(lp) {stats::dexp(x=x, rate=1/exp(lp))},
               numeric(length(x)))

      if(length(x) == 1) return(mean(pdmatrix))
      return(rowMeans(pdmatrix))
    }
  } else if(dist == "loglogistic") {
    function(x) {
      survmatrix <-
        vapply(lpvec,
               function(lp) {eha::dllogis(x=x, shape = 1/exp(scale.ln), scale = exp(lp))},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else {
    stop("Distribution not defined")
  }
}


