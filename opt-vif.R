#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(knockoff)
library(MASS)
library(parallel)

N <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])
myseed <- as.numeric(args[4])
mycores <- args[5]
set.seed(myseed)
cat("\nN = "); cat(N)
cat("\np = "); cat(p)
cat("\nnsim = "); cat(nsim)
cat("\nseed = "); cat(myseed)
cat("\ndesired cores = "); cat(mycores);
cat("\ndetectCores() = "); cat(detectCores())
mycores <- max(c(mycores, detectCores() - 1))
options(cores=mycores)
options(mc.cores=mycores)
cat("\n--using "); cat(mycores); cat(" cores\n")


k <- min(c(p/2, 10))
A <- 3.5
FDR <- 0.2
kcols <- sample(1:p, size = k, replace = FALSE)
BETA <- vector('numeric', length=p)
BETA[kcols] <- sample(c(A, -A), size = k, replace = TRUE)

fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
tpr <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))

rho <- 0.3
SigmaGen <- matrix(0, nrow=p, ncol=p)
SigmaGen[lower.tri(SigmaGen, diag = F)] <- rho
SigmaGen <- SigmaGen + t(SigmaGen) + diag(p)
# top eigenvalue of Sigma: 1 + (p-1) * a
# other eigenvalues of Sigma: 1 - a 

get_viffun <- function(Sigma_inv, tol=1e-5){
  f <- function(svec){
    Smat <- diag(svec)
    V <- 2 * Smat - Smat %*% Sigma_inv %*% Smat
    return(ifelse(min(eigen(V)$values) < tol, Inf, sum(diag(solve(V)))))
    #return(sum(diag(solve(V))))
  }
  return(f)
}

get_gradfun <- function(Sigma_inv){
  f <- function(svec){
    Smat <- diag(svec)
    p <- length(svec)
    V <- 2*Smat - Smat %*% Sigma_inv %*% Smat
    Vi <- solve(V)
    # derUS <- 2 * diag(p) - (Smat %*% Sigma_inv + Sigma_inv %*% Smat)
    sapply(1:p, function(j){
      Jjj <- matrix(0, nrow=p, ncol=p)
      Jjj[j,j] <- 1
      derUSjj <- 2 * Jjj - (Smat %*% Sigma_inv %*% Jjj + Jjj %*% Sigma_inv %*% Smat)
      sum(diag(t(-Vi %*% Vi) %*% derUSjj))
    })
  }
  return(f)
}


get_knockoffs <- function(Smat, Sigma_inv, X,Xsvd) {
  Cmat <- chol((N - 1 ) * (2 * Smat - Smat %*% Sigma_inv %*% Smat), pivot=T)
  pivot <- attr(Cmat,'pivot')
  po <- order(pivot)
  Cmat <- Cmat[,po]
  U<- Xsvd$u
  Q <- qr.Q(qr(cbind(U, matrix(0, nrow=N,ncol=p))))
  Utilde <- Q[,(p+1):(2*p)]
  ## Random Utilde
  Utilde <- Utilde %*% qr.Q(qr(matrix(rnorm(n=p*p), nrow=p, ncol=p)))
  return(X %*% (diag(p) - Sigma_inv %*% Smat) + Utilde %*% Cmat)
}


onesimrun <- function(){
  X <- mvrnorm(n=N, mu=rep(0, p), Sigma=SigmaGen)
  X <- scale(X, center=T) #scale(X, center=F, scale=sqrt(colSums(X^2)))
  Sigma <- cov(X)  #crossprod(X)
  Sigma_inv <- solve(Sigma)
  Xsvd <- svd(X)
  S_equi <- diag(rep(min(c(2 * min(eigen(Sigma)$values), 1))  , p))
  S_sdp <- diag(create.solve_sdp(Sigma))
  viffun <- get_viffun(Sigma_inv)
  
  opt_vif <- 
    optim(seq(0.1,0.9,length.out=p),
          fn=viffun, gr=get_gradfun(Sigma_inv),
          method = 'BFGS')
  #lower = rep(0, p),
  #upper = rep(1, p),
  #method='L-BFGS-B')
  S_vif <- diag(opt_vif$par)
  Y <- rnorm(N) + X %*% BETA
  Slist <- setNames( list(S_equi, S_sdp, S_vif),
                     c('equi', 'sdp', 'vif'))
  Xtlist <- lapply(Slist, function(SS) return(get_knockoffs(SS, Sigma_inv, X, Xsvd)))
  Xauglist <- lapply(Xtlist, function(Xtilde) return(cbind(X, Xtilde)))
  
  Wlist <- lapply(Slist, function(SS) return(vector('numeric', p)))
  blist <- lapply(Xauglist, function(Xaug) return(coef(lm(Y ~ 0 + Xaug))))
  Wlist <- lapply(blist, function(b) {
    W <- vector('numeric', p)
    for (j in seq_along(W)) {
      W[j] <- abs(b[j]) - abs(b[j + p])
    }
    return(W)
  })
  thresh <- lapply(Wlist, function(W) return(knockoff.threshold(W, fdr=FDR, offset=0)))
  sel <- mapply(W = Wlist, Tval = thresh,
                function(W, Tval) return(which(W >= Tval)),
                SIMPLIFY = FALSE)
  
  return(list(
    fdp = sapply(sel, function(ss) return(fdp(ss))),
    tpr = sapply(sel, function(ss) return(tpr(ss))),
    nsel = sapply(sel, length),
    vifsum =sapply(Slist, function(ss) return(viffun(diag(ss))))
  ))
}

res <- mclapply(1:nsim, function(i) return(onesimrun()))
save(res, N, p, BETA, k, nsim, myseed, mycores, A, FDR, rho, SigmaGen,
     file = paste('opt-vif-N', N, '-p', p, '.RData', sep=''))