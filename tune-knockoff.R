library(MASS)
get_viffun <- function(Sigma_inv, tol=1e-5){
  f <- function(svec){
    Smat <- diag(svec)
    V <- 2 * Smat - Smat %*% Sigma_inv %*% Smat
    return(ifelse(min(eigen(V)$values) < tol, Inf, sum(diag(solve(V)))))
    #return(sum(diag(solve(V))))
  }
  return(f)
}


get_all_svec <- function(Sigma){
  require(knockoff)
  list(
    sdp = knockoff::create.solve_sdp(Sigma),
    ldet = get_s_ldet(Sigma),
    equi = rep(min(c(1, 2 * min(eigen(Sigma, only.values=TRUE,symmetric = TRUE)$values))), ncol(Sigma))
  )
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

get_s_ldet <- function(Sigma){
  p <- ncol(Sigma)
  ldetGfunc <- get_ldetfun(Sigma)
  ldetGgrad <- get_ldetgrad(Sigma)
  s_init <- constrOptim(rep(0.001, p),
                        f = ldetGfunc,
                        grad = NULL,
                        ui = rbind(diag(p), -diag(p)),
                        ci = c(rep(0, p), rep(-1, p)),
                        control = list(fnscale = -1, maxit=5))$par
  Gopt <-
    constrOptim(s_init,
                f = ldetGfunc,
                grad = ldetGgrad,
                ui = rbind(diag(p), -diag(p)),
                ci = c(rep(0,p), rep(-1, p)),
                control = list(fnscale = -1))
  return(Gopt$par)
}

get_Cmat_eigen <- function(X, svec, G = NULL, Ginv = NULL, tol=1e-7){
  if (is.null(G)) {
    G <- crossprod(X)
  }
  if (is.null(Ginv)){
    Ginv <- solve(G)
  }
  Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
  CtC <- 2 * diag(svec) - diag(svec) %*% Ginv_S
  CtCeig <- eigen(CtC, symmetric=T)
  CtCeig$values[CtCeig$values < tol] <- 0
  Cmat <- diag(sqrt(CtCeig$values)) %*% t(CtCeig$vectors)
  return(Cmat)
}

get_knockoffs_qr <- function(X, svec, xqr = NULL, random = TRUE, tol = 1e-7, Cmat=NULL, Ginv = NULL, G= NULL){
  N <- nrow(X)
  p <- ncol(X)
  if (is.null(G)){
    G <- crossprod(X)
  }
  if (is.null(Cmat)){
    Cmat <- get_Cmat_eigen(X, svec, G, tol=tol)
  }
  if (is.null(Ginv)){
    Ginv <- solve(crossprod(X))
  }
  Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
  if (random){
    if (is.null(xqr)){
      Q <- qr.Q(qr(X))
    } else {
      Q <- qr.Q(xqr)
    }
    Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde <- Utilde - Q %*% crossprod(Q, Utilde) #(I - QQ^t) Utilde
    Utilde <- qr.Q(qr(Utilde))
  } else{
    Q <- qr.Q(qr(cbind(X, matrix(0, nrow=N, ncol=p))))
    Utilde <- Q[,(p+1):(2*p)]
  }
  return(X - X %*% Ginv_S + Utilde %*% Cmat)
}

get_knockoffs <- function(svec, X, Xsvd, random=FALSE) {
  N <- nrow(X)
  p <- ncol(X)
  VDi <- sweep(Xsvd$v, MARGIN=2, 1/Xsvd$d, `*`)
  SVDi <- sweep(VDi, MARGIN=1, svec, `*`)
  Sigma_inv_S <- tcrossprod(VDi, SVDi)
  nSSigiS <- -tcrossprod(SVDi)
  diag(nSSigiS) <- diag(nSSigiS) + 2 * svec
  Cmat <- chol(nSSigiS, pivot=T)
  pivot <- attr(Cmat,'pivot')
  po <- order(pivot)
  Cmat <- Cmat[,po]
  if (!random){
    Q <- qr.Q(qr(cbind(Xsvd$u, matrix(0, nrow=N,ncol=p))))
    Utilde <- Q[,(p+1):(2*p)]
  } else{
    Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde <- Utilde - Xsvd$u %*% crossprod(Xsvd$u, Utilde) #(I - UU^t) Utilde
    Utilde <- qr.Q(qr(Utilde))
  }
  return(X - X %*% Sigma_inv_S + Utilde %*% Cmat)
}


stat.ols.ginv <- function(X, Xk, y){
  XXk <- cbind(X,Xk) # N times 2*p
  b <- .lm.fit(XXk, y)$coefficients
  if (any(is.na(b))){
    XXk_gi <- ginv(XXk) # 2*p times N
    b <- XXk_gi %*% y
  }
  p <- ncol(X)
  W <- vector('numeric', p)
  for (j in seq_along(W)) {
    W[j] <- abs(b[j]) - abs(b[j + p])
  }
  return(W)
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
