library(knockoff)
library(MASS)


get_detfun <- function(Sigma){
  f <- function(svec){
    return(
      prod(svec) * det(2 * Sigma - diag(svec))
    )
  }
  return(f)
}

get_detgrad <- function(Sigma){
  p <- ncol(Sigma)
  f <- function(svec){
    S <- diag(svec)
    W <- 2 * Sigma - S
    Wdet <- det(W)
    ret <- vector('numeric', p)
    for (j in 1:p){
      Jjj <- matrix(0,nrow=p,ncol=p)
      Jjj[j,j] <- 1
      detWgrad <- Wdet * sum(diag(solve(W, -Jjj)))
      ret[j] <- prod(svec) * detWgrad + Wdet * prod(svec[-j])
    }
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

stat.olsdiff <- function(X, Xk, y){
  b <- .lm.fit(cbind(X,Xk), y)$coefficients
  p <- ncol(X)
  W <- vector('numeric', p)
  for (j in seq_along(W)) {
    W[j] <- abs(b[j]) - abs(b[j + p])
  }
  return(W)
}

canonical_svd <- function (X) {
  X.svd = tryCatch({
    svd(X)
  }, warning = function(w) {
  }, error = function(e) {
    stop("SVD failed in the creation of fixed-design knockoffs. Try upgrading R to version >= 3.3.0")
  }, finally = {
  })
  for (j in 1:min(dim(X))) {
    i = which.max(abs(X.svd$u[, j]))
    if (X.svd$u[i, j] < 0) {
      X.svd$u[, j] = -X.svd$u[, j]
      X.svd$v[, j] = -X.svd$v[, j]
    }
  }
  return(X.svd)
}

normc <- function (X, center = T) {
  X.centered = scale(X, center = center, scale = F)
  X.scaled = scale(X.centered, center = F, scale = sqrt(colSums(X.centered^2)))
  X.scaled[, ]
}

get_Xgenfunc <- function(SigmaGen, N){
  L <- chol(SigmaGen)
  p <- ncol(SigmaGen)
  f <- function() {
      matrix(rnorm(N*p), nrow=N, ncol=p) %*% L
  }
  return(f)
}

onesimrun <- function(SigmaGen, Xgenfunc, BETA, N, FDR, statfunclist) {
  statnames <- names(statfunclist)
  p <- ncol(SigmaGen)
  X <- Xgenfunc()
  X <- normc(X, center=F)
  Sigma <- crossprod(X)
  Xsvd <- canonical_svd(X)
  svec_equi <- rep(min(c(2 * min(Xsvd$d^2), 1))  , p)
  if (p > 100){
    svec_sdp <- create.solve_asdp(Sigma)
  } else{
    svec_sdp <- create.solve_sdp(Sigma)
  }
  ldetGfunc <- get_ldetfun(Sigma)
  ldetGgrad <- get_ldetgrad(Sigma)
  # Optimize over 0 <= s_j <= 1
  Gopt <-
    constrOptim(rep(0.001,p),
                f = ldetGfunc,
                grad = ldetGgrad,
                ui = rbind(diag(p), -diag(p)),
                ci = c(rep(0,p), rep(-1, p)),
                control = list(fnscale = -1))
  Y <- rnorm(N) + X %*% BETA
  sveclist <- setNames(
    list(svec_equi, svec_sdp, Gopt$par),
    c('equi', 'sdp', 'Gdet')
  )
  Xtlist <- lapply(sveclist, function(ss) return(get_knockoffs(ss, X, Xsvd)))
  Wlist <- 
    lapply(statfunclist, function(stf) return(
      lapply(Xtlist, function(Xk) return(stf(X, Xk, Y)))
  ))
  thresh <- 
    lapply(Wlist,
         function(statW) return(
           lapply(statW, function(W) return(knockoff.threshold(W, fdr=FDR, offset=1)))
         ))
  sel <- 
    mapply(W_stat = Wlist, thresh_stat = thresh,
           function(W_stat, thresh_stat){
             mapply(W = W_stat, Tval = thresh_stat,
                    function(W, Tval) return(which(W >= Tval)),
                    SIMPLIFY = FALSE)
           }, SIMPLIFY = F
           )
  ret <- lapply(sel,
         function(sel_stat) return(
           list(
             fdp = sapply(sel_stat, function(ss) return(fdp(ss))),
             tpr = sapply(sel_stat, function(ss) return(tpr(ss))),
             ppv = sapply(sel_stat, function(ss) return(ppv(ss))),
             nsel = sapply(sel_stat, length),
             lGdet = sapply(sveclist, function(ss) return(ldetGfunc(ss))),
             s = do.call('cbind', sveclist)
           )
         ))
  for(j in seq_along(ret)){
    ret[[j]]$statname <- statnames[j]
  }
  return(ret)
}