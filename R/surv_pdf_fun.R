
# Create a function to calculate S(t)
# The function to be created will take time (`x`) as the input argument and return
create_survfun <- function(lpvec, scale,
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
               function(lp) {stats::plnorm(q=x, meanlog=lp, sdlog=exp(scale), lower=FALSE)},
               numeric(length(x)))

      if(length(x) == 1) return(mean(survmatrix))
      return(rowMeans(survmatrix))
    }
  } else {
    stop("Distribution not defined")
  }
}

# Create a function to calculate PDF(t)
create_pdf <- function(lpvec, scale, dist = "lognormal"){
  function(x) {
    pdmatrix <-
      vapply(lpvec,
             function(lp) {stats::dlnorm(x=x, meanlog=lp, sdlog=exp(scale))},
             numeric(length(x)))

    if(length(x) == 1) return(mean(pdmatrix))
    return(rowMeans(pdmatrix))
  }
}


