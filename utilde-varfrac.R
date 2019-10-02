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
rho <- as.numeric(args[10])
statname <- args[11]
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
if (statname == 'cordiff'){
  cat("\nW statistic: cross product differences")
  wstatfunc <- stat.crossprod
} else if (statname == 'lasso-lambda'){
  wstatfunc <- stat.lasso_lambdadiff
  cat("\nW statistic: difference in lasso lambda")
} else if (statname == 'ridge'){
  wstatfunc <- stat.ridge
  cat("\nW statistic: ridge regression")
} else {
  cat("\nunknown W statistic, using cross products")
  wstatfunc <- stat.crossprod
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
simres <- 
  foreach(i = 1:nsim, .combine = rbind) %dopar% {
    ret_i <- matrix(nrow = 2, ncol = result_length)
    colnames(ret_i) <- resnames
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
    Utilde <- get_Utilde_random(Qx, N, p)
    Utheta <- get_Utheta(Qx, Y, Yresid_normed, N, p, ufrac * ufrac_multiplier)
    Xtilde_Utheta <-  X_minus_XGinvS + Utheta %*% Cmat
    Xtilde <- X_minus_XGinvS + Utilde %*% Cmat
    W_Utheta <- wstatfunc(X, Xtilde_Utheta, Y)
    W_regular <- wstatfunc(X, Xtilde, Y)
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


save(res_fmt, BETA, myseed, file = paste("utheta-gridsearch-stat", statname, "-N", N, "-p", p, "-rho", rho,
                                         "-off", offset, '-',sigmatype,
                                         '.RData', 
                                         sep='')
)