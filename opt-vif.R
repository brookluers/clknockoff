get_viffun <- function(Sigma_inv, tol=1e-5){
  f <- function(svec){
    Smat <- diag(svec)
    V <- 2 * Smat - Smat %*% Sigma_inv %*% Smat
    return(ifelse(min(eigen(V)$values) < tol, Inf, sum(diag(solve(V)))))
    #return(sum(diag(solve(V))))
  }
  return(f)
}

get_vifgrad <- function(Sigma_inv){
  p <- ncol(Sigma_inv)
  Jlist <- vector('list', p)
  for (k in seq_along(Jlist)){
    Jlist[[k]] <- matrix(0,nrow=p,ncol=p)
    Jlist[[k]][k,k] <- 1
  }
  f <- function(svec){
    Smat <- diag(svec)
    U <- 2 * Smat - Smat %*% Sigma_inv %*% Smat
    Ui <- solve(U)
    ret <- sapply(1:p, function(j){
      dUds <- matrix(0, nrow=p, ncol=p)
      dUds[,j] <- - ( Sigma_inv[,j] * svec )
      dUds[j,] <- dUds[j,] - svec * Sigma_inv[j,]
      dUds[j,j] <- dUds[j,j] + 2
      return(sum(diag( -Ui %*% dUds %*% Ui )))
    })
    return(ret)
  }
  return(f)
}