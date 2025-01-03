int2BinaryStr <- function(x) {
  min <- -2^31 # min value for 32 bit integer
  max <- 2^31 - 1 # max value for 32 bit integer
  if (is.integer(x) && x >= min && x <= max){
    bits <- intToBits(as.integer(x))
    bin_str <- paste(rev(as.integer(bits)), collapse = "")

    return(bin_str)
  }
  else {
    return(NA)
  }
}