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

get_ldetfun <- function(Sigma) {
  # det(G) = det(S) * det(2Sigma - S)
  #log(det(G)) = sum(log(eigenvalues of S)) + sum(log(eigenvalues of 2Sigma - S))
  f <- function(svec) {
    # Wdet <- determinant(2 * Sigma - diag(svec), logarithm = T)
    W <- 2 * Sigma - diag(svec)
    Wev <- eigen(W, symmetric = T, only.values = T)$values
    if (any(Wev < 0)) {
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

get_knockoffs <- function(Smat, Sigma_inv, X, Xsvd) {
  Cmat <- chol(2 * Smat - Smat %*% Sigma_inv %*% Smat, pivot=T)
  pivot <- attr(Cmat,'pivot')
  po <- order(pivot)
  Cmat <- Cmat[,po]
  U <- Xsvd$u
  Q <- qr.Q(qr(cbind(U, matrix(0, nrow=N,ncol=p))))
  Utilde <- Q[,(p+1):(2*p)]
  ## Random Utilde
  #Utilde <- Utilde %*% qr.Q(qr(matrix(rnorm(n=p*p), nrow=p, ncol=p)))
  return(X %*% (diag(p) - Sigma_inv %*% Smat) + Utilde %*% Cmat)
}

stat.olsdiff <- function(X, Xk, y){
  b <- coef(lm(y ~ 0 + cbind(X, Xk)))
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

onesimrun <- function(SigmaGen, BETA, N, FDR, statfunc = stat.olsdiff, statname='ols'){
  p <- ncol(SigmaGen)
  X <- mvrnorm(n=N, mu=rep(0, p), Sigma=SigmaGen)
  # X <- scale(X, center=T) 
  X <- normc(X, center=F)
  #Sigma <- cov(X)
  Sigma <- crossprod(X)
  Sigma_inv <- solve(Sigma)
  Xsvd <- canonical_svd(X)
  S_equi <- diag(rep(min(c(2 * min(eigen(Sigma)$values), 1))  , p))
  S_sdp <- diag(create.solve_sdp(Sigma))
  
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
  S_Gdet <- diag(Gopt$par)
  Y <- rnorm(N) + X %*% BETA
  Slist <- setNames( list(S_equi, S_sdp, S_Gdet),
                     c('equi', 'sdp', 'Gdet'))
  Xtlist <- lapply(Slist, function(SS) return(get_knockoffs(SS, Sigma_inv, X, Xsvd)))
  Xauglist <- lapply(Xtlist, function(Xtilde) return(cbind(X, Xtilde)))
  Glist <- lapply(Xauglist, crossprod)
  Geig <- lapply(Glist, function(G) return(eigen(G,only.values=T,symmetric=T)$values))
  Wlist <- lapply(Xtlist, function(Xk) return(statfunc(X, Xk, Y)))
  thresh <- lapply(Wlist, function(W) return(knockoff.threshold(W, fdr=FDR, offset=1)))
  sel <- mapply(W = Wlist, Tval = thresh,
                function(W, Tval) return(which(W >= Tval)),
                SIMPLIFY = FALSE)
  
  return(list(
    fdp = sapply(sel, function(ss) return(fdp(ss))),
    tpr = sapply(sel, function(ss) return(tpr(ss))),
    ppv = sapply(sel, function(ss) return(ppv(ss))),
    nsel = sapply(sel, length),
    Geig = sapply(Geig,function(x)return(c(min(x),max(x)))),
    lGdet = sapply(Slist, function(ss) return(ldetGfunc(diag(ss)))),
    s = do.call('cbind',lapply(Slist, function(ss)return(diag(ss)))),
    statname = statname
  ))
}