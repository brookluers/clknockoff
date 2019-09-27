#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(doParallel)
library(knockoff)
library(MASS)
source("tune-knockoff.R")
source("pop-cov-funcs.R")
N <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])
myseed <- as.numeric(args[4])
mycores <- as.numeric(args[5])
sigmatype <- args[6]
k <- as.numeric(args[7])
offset <- as.numeric(args[8])
ufrac_multiplier <- as.numeric(args[9])
if (is.na(args[10])){
  rho <- 0.8
} else {
  rho <- as.numeric(args[10])
}
RNGkind("L'Ecuyer-CMRG")
set.seed(myseed)
cat("\nN = "); cat(N)
cat("\np = "); cat(p)
cat("\nnsim = "); cat(nsim)
cat("\nseed = "); cat(myseed)
cat("\ndesired cores = "); cat(mycores);
cat("\ndetectCores() = "); cat(detectCores())
mycores <- max(c(mycores - 1, detectCores() - 1))
options(cores=mycores)
registerDoParallel(cores = mycores)
cat("\n--using "); cat(mycores); cat(" cores\n")
cat("\nMultiply target ||UU^t Y|| fraction by "); cat(ufrac_multiplier); cat("\n")
if (sigmatype == 'exch'){
  SigmaGen <- get_exch(p, rho)
} else if (sigmatype=='ar1'){
  SigmaGen <- get_ar(p, rho)
} else {
  cat("unknown Sigma structure\nusing exchangeable\n\n")
  SigmaGen <- get_exch(p, rho)
}

norm_utheta_projy <- function(theta, Utilde1, Utilde2, Y){
  ret <- vector('numeric', length(theta))
  for (i in seq_along(ret)){
    Utheta <- sin(theta[i]) * Utilde1 + cos(theta[i]) * Utilde2
    ret[i] <- norm(Utheta %*% crossprod(Utheta, Y), 'F')^2
  }
  return(ret)
}
tseq <- seq((6/16)*pi, (10/16)*pi, length.out=200)


BETA <- vector('numeric', p)
fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
ppv <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))
tpr <- function(selected) sum(BETA[selected] != 0) / k
fpr <- function(selected) sum(BETA[selected] == 0) / (p - k)
FDR <- 0.1
magnitude <- 3.5
BETA[sample(1:p, size = k)] <- magnitude * sample(c(1,-1), size=k, replace=T)
one_to_p <- 1:p
selnames <- paste('sel', 1:p, sep='')
resnames <- c('sim_ix', 'uuynorm', 'uuynormfrac', 'ufrac', 
              'fdp','fpr','ppv','tpr','nsel', selnames)
result_length <- length(resnames)
### Single X realization
X <- mvrnorm(N, mu = rep(0, p), Sigma=SigmaGen)
X <- scale(X, center = T, scale = F)
X <- scale(X, center = F, scale = apply(X, 2, function(xj)
  return(sqrt(sum(xj ^ 2)))))
Xbeta <- X %*% BETA
xqr <- qr(X)
Qx <- qr.Q(xqr)
G <- crossprod(X)
Ginv <- solve(G)
svec <- get_s_ldet(G)
Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
X_minus_XGinvS <- X - X %*% Ginv_S
Cmat <- get_Cmat_eigen(X, svec, G, Ginv)
### 
simres <- 
  foreach(i = 1:nsim, .combine = rbind) %dopar% {
    ret_i <- matrix(nrow = 2, ncol = result_length)
    colnames(ret_i) <- resnames
    Y <- Xbeta + rnorm(N)
    Ynorm2 <- norm(Y, 'F')^2
    XYcp <- crossprod(X, Y)
    abs_XYcp <- abs(XYcp)
    Ylm <- lm(Y ~ 0 + X)
    betahat <- coef(Ylm)
    Yresid <- residuals(Ylm)
    sig2hat <- sum(Yresid^2) / (N - p)
    ufrac <- (sig2hat * p) / ((N - p) * sig2hat + as.numeric(t(betahat) %*% G %*% cbind(betahat)))
    norm_Yresid <- norm(cbind(Yresid), type='F')
    Yresid_normed <- Yresid / norm_Yresid
    Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde, project away from X
    Utilde <- qr.Q(qr(Utilde)) # orthogonalize
    Q_xuy <- qr.Q(qr(cbind(Utilde, Yresid_normed, Qx)))
    Utilde_yperp <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde_yperp <- Utilde_yperp - Q_xuy %*% crossprod(Q_xuy, Utilde_yperp) # project away from Yresid, X, Utilde
    Utilde_yperp <- qr.Q(qr(Utilde_yperp)) # orthogonalize
    Utilde_contain_yperp <- cbind(Yresid_normed, matrix(rnorm(N*(p-1)), nrow=N, ncol=(p-1)))
    Q_xu <- qr.Q(qr(cbind(Utilde, Qx)))
    Utilde_contain_yperp <- Utilde_contain_yperp - Q_xu %*% crossprod(Q_xu, Utilde_contain_yperp)
    Utilde_contain_yperp <- qr.Q(qr(Utilde_contain_yperp))
    nu2y_tseq <- norm_utheta_projy(tseq, Utilde, Utilde_yperp, Y)
    nu2yfrac_tseq <- nu2y_tseq / Ynorm2
    nu3y_tseq <- norm_utheta_projy(tseq, Utilde, Utilde_contain_yperp, Y)
    nu3yfrac_tseq <- nu3y_tseq / Ynorm2
    minu3frac.ix <- which.min(abs(nu3yfrac_tseq - ufrac_multiplier * ufrac))
    minu2frac.ix <- which.min(abs(nu2yfrac_tseq - ufrac_multiplier * ufrac))
    if (abs(nu3yfrac_tseq[minu3frac.ix] - ufrac_multiplier * ufrac) < abs(nu2yfrac_tseq[minu2frac.ix] - ufrac_multiplier * ufrac)) {
      # Use Utilde_3 definition, = Utilde_contain_yperp
      theta <- tseq[minu3frac.ix]
      Utheta <- sin(theta) * Utilde + cos(theta) * Utilde_contain_yperp
    } else {
      # Use Utilde_2   = Utilde_yperp
      theta <- tseq[minu2frac.ix]
      Utheta <- sin(theta) * Utilde + cos(theta) * Utilde_yperp
    }
    
    Xtilde_Utheta <-  X_minus_XGinvS + Utheta %*% Cmat
    Xtilde <- X_minus_XGinvS + Utilde %*% Cmat
    W_Utheta <- as.numeric(abs_XYcp - abs(crossprod(Xtilde_Utheta, Y))) # simple correlation differences
    W_regular <- as.numeric(abs_XYcp - abs(crossprod(Xtilde, Y)))
    sel_Utheta <- which(W_Utheta >= knockoff.threshold(W_Utheta, FDR, offset=offset))
    sel_regular <- which(W_regular >= knockoff.threshold(W_regular, FDR, offset=offset))
    
    ret_i <- 
      rbind(
        c(i, norm(Utheta %*% crossprod(Utheta, Y), 'F')^2,
          NA, # uunormfrac
          ufrac, 
          fdp(sel_Utheta),
          fpr(sel_Utheta),  ppv(sel_Utheta),
          tpr(sel_Utheta), length(sel_Utheta),
          1 * (one_to_p %in% sel_Utheta)
        ),
        c(i, norm(Utilde %*% crossprod(Utilde, Y), 'F')^2,
          NA, # uunormfrac
          ufrac, 
          fdp(sel_regular),
          fpr(sel_regular),  ppv(sel_regular),
          tpr(sel_regular), length(sel_regular),
          1 * (one_to_p %in% sel_regular)
        )
      )
    colnames(ret_i) <- resnames
    rownames(ret_i) <- c('Utheta', 'Utilde')
    ret_i[,'uuynormfrac'] <- ret_i[,'uuynorm'] / Ynorm2
    ret_i
  }

res_fmt <- as_tibble(simres, rownames = 'utype')
res_fmt$pop_Sigma_cnum <- kappa(SigmaGen, exact = TRUE)
res_fmt$k <- k
res_fmt$FDR <- FDR
res_fmt$signal <- magnitude
res_fmt$sigmatype <- sigmatype
res_fmt$rho <- rho
res_fmt$offset <- offset
res_fmt$ufrac_multiplier <- ufrac_multiplier
res_fmt$N <- N
res_fmt$p <- p


save(res_fmt, X, BETA, myseed, file = paste("utheta-fixX-mult", ufrac_multiplier, "-N", N, "-p", 
                                         p, "-rho", rho,
                                         "-seed", myseed,
                                         "-off", offset, '-', sigmatype,
                                         '.RData', 
                                         sep='')
)