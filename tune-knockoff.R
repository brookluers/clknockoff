
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

get_ldetfun <- function(Sigma, tol=1e-12) {
  # det(G) = det(S) * det(2Sigma - S)
  #log(det(G)) = sum(log(eigenvalues of S)) + sum(log(eigenvalues of 2Sigma - S))
  f <- function(svec) {
    # Wdet <- determinant(2 * Sigma - diag(svec), logarithm = T)
    W <- 2 * Sigma - diag(svec)
    Wev <- eigen(W, symmetric = T, only.values = T)$values
    if (any(Wev < tol)) {
      return(-Inf)
    } else{
      return(
        # sum(log(svec)) + Wdet$modulus[1]
        sum(log(svec)) + sum(log(Wev))
      )
    }
    
  }
  return(f)
}

get_ldetgrad <- function(Sigma) {
  p <- ncol(Sigma)
  ## single-entry matrix (jj) selects the jth column
  # when right multiplied
  f <- function(svec){
    ret <- vector('numeric', p)
    W <- 2 * Sigma - diag(svec)
    Winv <- solve(W)
    ret <- 1 / svec - diag(Winv)
    return(ret)
  }
  return(f)
}


get_knockoffs_old <- function(Smat, Sigma_inv, X, Xsvd) {
  Cmat <- chol(2 * Smat - Smat %*% Sigma_inv %*% Smat, pivot=T)
  pivot <- attr(Cmat,'pivot')
  po <- order(pivot)
  Cmat <- Cmat[,po]
  U <- Xsvd$u
  Q <- qr.Q(qr(cbind(U, matrix(0, nrow=N,ncol=p))))
  Utilde <- Q[,(p+1):(2*p)]
  ## Random Utilde
  #  Utilde <- Utilde %*% qr.Q(qr(matrix(rnorm(n=p*p), nrow=p, ncol=p)))
  return(X %*% (diag(p) - Sigma_inv %*% Smat) + Utilde %*% Cmat)
}

get_knockoffs <- function(svec, X, Xsvd) {
  VDi <- sweep(Xsvd$v, MARGIN=2, 1/Xsvd$d, `*`)
  SVDi <- sweep(VDi, MARGIN=1, svec, `*`)
  Sigma_inv_S <- tcrossprod(VDi, SVDi)
  nSSigiS <- -tcrossprod(SVDi)
  diag(nSSigiS) <- diag(nSSigiS) + 2 * svec
  Cmat <- chol(nSSigiS, pivot=T)
  pivot <- attr(Cmat,'pivot')
  po <- order(pivot)
  Cmat <- Cmat[,po]
  U <- Xsvd$u
  Q <- qr.Q(qr(cbind(U, matrix(0, nrow=N,ncol=p))))
  Utilde <- Q[,(p+1):(2*p)]
  return(X - X %*% Sigma_inv_S + Utilde %*% Cmat)
}

stat.olsdiff <- function(X, Xk, y){
  b <- .lm.fit(cbind(X,Xk), y)$coefficients
  p <- ncol(X)
  W <- vector('numeric', p)
  for (j in seq_along(W)) {
    W[j] <- abs(b[j]) - abs(b[j + p])
  }
  return(W)
}

stat.ridge <- function(X, X_k, y, lambda=NULL){
  if (is.null(lambda)){
    lambda <- c(seq(0.01,0.99,length.out=20))
  }
  p <- ncol(X)
  rfit <- lm.ridge(y ~ 0 + cbind(X, X_k), lambda=lambda)
  if( length(lambda)>1){
    b <- coef(rfit)[which.min(rfit$GCV),]
  } else {
    b <- coef(rfit)
  }
  
  W <- vector('numeric', p)
  for (j in seq_along(W)){
    W[j] <- abs(b[j]) - abs(b[j + p])
  }
  return(W)
}