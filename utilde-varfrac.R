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
if (is.na(args[9])){
  rho <- 0.8
} else {
  rho <- as.numeric(args[9])
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

if (sigmatype == 'exch'){
  SigmaGen <- get_exch(p, rho)
} else if (sigmatype=='ar1'){
  SigmaGen <- get_ar(p, rho)
} else {
  cat("unknown Sigma structure\nusing exchangeable\n\n")
  SigmaGen <- get_exch(p, rho)
}

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
resnames <- c('sim_ix', 'uuynorm', 'uuynormfrac', 'fdp','fpr','ppv','tpr','nsel', selnames)
result_length <- length(resnames)
simres <- 
  foreach(i = 1:nsim, .combine = rbind) %dopar% {
    ret_i <- matrix(nrow = 2, ncol = result_length)
    colnames(ret_i) <- resnames
    rownames(ret_i) <- c("Uavg", "Utilde")
    X <- mvrnorm(N, mu = rep(0,p),Sigma=SigmaGen)
    X <- scale(X, center = T, scale = F)
    X <- scale(X, center = F, scale = apply(X, 2, function(xj)
      return(sqrt(sum(xj ^ 2)))))
    Xbeta <- X %*% BETA
    offset <- 1
    xqr <- qr(X)
    Qx <- qr.Q(xqr)
    G <- crossprod(X)
    Ginv <- solve(G)
    svec <- get_s_ldet(G)
    Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
    X_minus_XGinvS <- X - X %*% Ginv_S
    Cmat <- get_Cmat_eigen(X, svec, G, Ginv)
    Y <- Xbeta + rnorm(N)
    XYcp <- crossprod(X, Y)
    abs_XYcp <- abs(XYcp)
    Ylm <- lm(Y ~ 0 + X)
    betahat <- coef(Ylm)
    ufrac <- p / (N + as.numeric(t(betahat) %*% G %*% cbind(betahat)))
    Yresid <- residuals(Ylm)
    norm_Yresid <- norm(cbind(Yresid), type='F')
    Yresid_normed <- Yresid / norm_Yresid
    Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde, project away from X
    Utilde <- qr.Q(qr(Utilde)) # orthogonalize
    Q_xuy <- qr.Q(qr(cbind(Utilde, Yresid_normed, Qx)))
    Utilde_yperp <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde_yperp <- Utilde_yperp - Q_xuy %*% crossprod(Q_xuy, Utilde_yperp) # project away from Yresid, X, Utilde
    Utilde_yperp <- qr.Q(qr(Utilde_yperp)) # orthogonalize
    Utilde_avg  <- ufrac * Utilde + (1 - ufrac) * Utilde_yperp
    Xtilde_Uavg <-  X_minus_XGinvS + Utilde_avg %*% Cmat
    Xtilde <- X_minus_XGinvS + Utilde %*% Cmat
    W_Uavg <- as.numeric(abs_XYcp - abs(crossprod(Xtilde_Uavg, Y))) # simple correlation differences
    W_regular <- as.numeric(abs_XYcp - abs(crossprod(Xtilde, Y)))
    sel_Uavg <- which(W_Uavg >= knockoff.threshold(W_Uavg, FDR, offset=offset))
    sel_regular <- which(W_regular >= knockoff.threshold(W_regular, FDR, offset=offset))
    
    ret_i <- 
      rbind(
        c(i, norm(Utilde_avg %*% crossprod(Utilde_avg, Y), 'F')^2,
          NA, # uunormfrac
          fdp(sel_Uavg),
          fpr(sel_Uavg),  ppv(sel_Uavg),
          tpr(sel_Uavg), length(sel_Uavg),
          1 * (one_to_p %in% sel_Uavg)
        ),
        c(i, norm(Utilde %*% crossprod(Utilde, Y), 'F')^2,
          NA, # uunormfrac
          fdp(sel_regular),
          fpr(sel_regular),  ppv(sel_regular),
          tpr(sel_regular), length(sel_regular),
          1 * (one_to_p %in% sel_regular)
        )
      )
    colnames(ret_i) <- resnames
    rownames(ret_i) <- c('Uavg', 'Utilde')
    ret_i[,'uuynormfrac'] <- ret_i[,'uuynorm'] / norm(Y, 'F')^2
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
res_fmt$N <- N
res_fmt$p <- p


save(res_fmt, BETA, myseed, file = paste("utilde-fraction-statcordiff-N", N, "-p", p, "-rho", rho,
                                         "-off", offset, '-',sigmatype,
                                         '.RData', 
                                         sep='')
)