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
    return(1 / svec - diag(solve(2 * Sigma - diag(svec))))
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

get_Utilde_random <- function(Qx, N, p){
  Utilde <- matrix(rnorm(N * p), nrow = N, ncol = p)
  Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde
  return(qr.Q(qr(Utilde)))
}

norm_utheta_projy <- function(theta, Utilde1, Utilde2, Y){
  ret <- vector('numeric', length(theta))
  for (i in seq_along(ret)){
    Utheta <- sin(theta[i]) * Utilde1 + cos(theta[i]) * Utilde2
    ret[i] <- norm(Utheta %*% crossprod(Utheta, Y), 'F')^2
  }
  return(ret)
}


get_Utheta <-
  function(Qx, Y, Yresid_normed,  N, p, target_val, Ynorm2 = NULL, tseq = NULL) {
    if (is.null(tseq)) {
      tseq <- seq((1 / 4) * pi, (3 / 4) * pi, length.out = 250)
    }
    if (is.null(Ynorm2)){
      Ynorm2 <-  norm(Y, 'F')^2
    }
    Utilde <- get_Utilde_random(Qx, N, p)
    Q_xuy <- qr.Q(qr(cbind(Utilde, Yresid_normed, Qx)))
    Utilde_yperp <- matrix(rnorm(N * p), nrow = N, ncol = p)
    Utilde_yperp <-
      Utilde_yperp - Q_xuy %*% crossprod(Q_xuy, Utilde_yperp) # project away from Yresid, X, Utilde
    Utilde_yperp <- qr.Q(qr(Utilde_yperp)) # orthogonalize
    Utilde_contain_yperp <-
      cbind(Yresid_normed, matrix(rnorm(N * (p - 1)), nrow = N, ncol = (p - 1)))
    Q_xu <- qr.Q(qr(cbind(Utilde, Qx)))
    Utilde_contain_yperp <-
      Utilde_contain_yperp - Q_xu %*% crossprod(Q_xu, Utilde_contain_yperp)
    Utilde_contain_yperp <- qr.Q(qr(Utilde_contain_yperp))
    
    nu2y_tseq <- norm_utheta_projy(tseq, Utilde, Utilde_yperp, Y)
    nu2yfrac_tseq <- nu2y_tseq / Ynorm2
    nu3y_tseq <-
      norm_utheta_projy(tseq, Utilde, Utilde_contain_yperp, Y)
    nu3yfrac_tseq <- nu3y_tseq / Ynorm2
    minu3frac.ix <-
      which.min(abs(nu3yfrac_tseq - target_val))
    minu2frac.ix <-
      which.min(abs(nu2yfrac_tseq - target_val))
    if (abs(nu3yfrac_tseq[minu3frac.ix] - target_val) < abs(nu2yfrac_tseq[minu2frac.ix] - target_val)) {
      # Use Utilde_3 definition, = Utilde_contain_yperp
      theta <- tseq[minu3frac.ix]
      Utheta <-
        sin(theta) * Utilde + cos(theta) * Utilde_contain_yperp
    } else {
      # Use Utilde_2   = Utilde_yperp
      theta <- tseq[minu2frac.ix]
      Utheta <- sin(theta) * Utilde + cos(theta) * Utilde_yperp
    }
    return(Utheta)
  }

get_knockoffs_qr <- function(X, svec, xqr = NULL, utilde = 'random', tol = 1e-7, Cmat=NULL, Ginv = NULL, G= NULL,
                             ufrac_target = NULL, Y = NULL){
  N <- nrow(X)
  p <- ncol(X)
  if (is.null(G)){
    G <- crossprod(X)
  }
  if (is.null(Cmat)){
    Cmat <- get_Cmat_eigen(X, svec, G, tol=tol)
  }
  if (is.null(Ginv)){
    Ginv <- solve(G)
  }
  if (is.null(xqr)){
    xqr <- qr(X)
  }
  Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
  if (utilde == 'random'){
    Qx <- qr.Q(xqr)
    Utilde <- get_Utilde_random(Qx, N, p)
  } else if (utilde=='utheta') {
    Yresid <- residuals(lm(Y ~ 0 + X))
    norm_Yresid <- norm(cbind(Yresid), type='F')
    Yresid_normed <- Yresid / norm_Yresid
    Qx <- qr.Q(xqr)
    Utilde <- get_Utheta(Qx, Y, Yresid_normed,N,p,ufrac_target,Ynorm2 = norm(Y,'F')^2)
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

stat.crossprod <- function(X, Xk, y){
  XYcp <- crossprod(X, y)
  abs_XYcp <- abs(XYcp)
  return(as.numeric(abs_XYcp - abs(crossprod(Xk, y))))
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
