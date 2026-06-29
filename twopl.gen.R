twopl.gen <- function(a,b,theta){
  
  twopl <- function(a,b,theta){
    num <- exp(a*(theta-b))
    deno <- 1+exp(a*(theta-b))
    prob <- num/deno
    return(prob)
  }
  
  J <- length(a)
  N <- length(theta)
  dat <- matrix(NA, nr=N, nc=J)
  rand <- matrix(runif(N*J,0,1), nr=N, nc=J)
  
  for(i in 1:N){
    for(j in 1:J){
      tmp <- twopl(a[j], b[j], theta[i])
      if(tmp > rand[i,j]) dat[i,j] <- 1
      else dat[i,j] <- 0
    }
  }
  dat <- data.frame(dat)
  colnames(dat) <- paste0("Item",1:J)
  return(dat)
}