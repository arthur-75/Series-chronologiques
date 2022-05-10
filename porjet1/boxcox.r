boxcox <- function(x,a){
  if (a==0){
    res <- log(x)
  }
  else res <- (x^a-1)/a
  res <- ts(res,start=start(x),frequency=frequency(x))
}