get_exch <- function(p, rho){
  ret <- diag(p)
  ret[lower.tri(ret, diag = F)]<- rho
  ret[upper.tri(ret, diag=F)] <- rho
  return(ret)
}

get_ar <- function(p, rho){
  ret <- diag(p)
  ret <- rho^abs(row(ret) - col(ret))
  return(ret)
}

get_banded <- function(p, rho, nbands = 1){
  if (nbands > p-1) stop("nbands must be <= p - 1")
  ret <- diag(p)
  ixmat <- abs(row(ret) - col(ret))
  ret[ixmat <= 1 & ixmat > 0] <- rho
  return(ret)
}
