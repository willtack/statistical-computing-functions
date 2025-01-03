my_pexp <- function(x, rate, lower.tail=FALSE, log.p=FALSE) {
  
  if (is.nan(x) || is.nan(rate)){
    stop("At least one of x or rate is NaN.")
  }
  
  if (x < 0) {
    warning("x must be non-negative.")
    return(if (log.p) -Inf else 0)
  }
  
  if (rate <= 0) {
    stop("rate must be positive.")
  }
  
  pow <- -rate*x
  
  if (lower.tail){
    if(log.p){
      # handle ultra small numbers with approximation
      if (rate*x < 1e-10) {return(log(rate*x))} else {return(log1p(-exp(pow)))}
    } else {
      return(-expm1(pow))
    }
  } else {
    if(log.p){
      return(pow)
    } else {
      return(exp(pow))
    }
  }
}
