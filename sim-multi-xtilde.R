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
magnitude <- 3.5
BETA[sample(1:p, size = k)] <- magnitude * sample(c(1,-1), size=k, replace=T)
FDR <- 0.1
cat("Nominal FDR = "); cat(FDR); cat("\n")
fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
ppv <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))
tpr <- function(selected) sum(BETA[selected] != 0) / k
fpr <- function(selected) sum(BETA[selected] == 0) / (p - k)
one_to_p <- 1:p
selnames <- paste('sel', 1:p, sep='')
resnames <- c('fdp','fpr','ppv','tpr','nsel', selnames)
result_length <- length(resnames)
simres <- vector('list', nsim)
nxtil <- 50 # number of Xtilde matrices per fixed (X, Y)
pzeroes <- rep(0, p)
select_manyXtilde <- function(svec, X, Qx, Y, abs_XYcp, Cmat, Ginv, FDR, offset, nxtil = 10){
  Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
  X_minus_GinvS <- X - X %*% Ginv_S
  Wlist <- vector('list', nxtil)
  for(j in seq_along(Wlist)){
    Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde <- qr.Q(qr(Utilde)) # orthogonalize only, do not project away from proj_xperp(Y)
    Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde, project away from X
    Xtilde <-  X_minus_GinvS + Utilde %*% Cmat
    Wlist[[j]] <- as.numeric(abs_XYcp - abs(crossprod(Xtilde, Y)))
  }
  # choose vars s.t. estimated prob. of selecting them is at least FDR
  selmat <- do.call('rbind', lapply(Wlist, function(Wm) return(1 * (Wm >= knockoff.threshold(Wm, FDR, offset=offset)))))
  return(which(apply(selmat, 2, mean) >= FDR))
}

simres <- 
  foreach(i = 1:length(simres), .combine = rbind) %dopar% {
    X <- mvrnorm(N, mu = pzeroes, Sigma = SigmaGen)
    X <- scale(X, center = T, scale = F)
    X <-
      scale(X,
            center = F,
            scale = apply(X, 2, function(xj)
              return(sqrt(sum(
                xj ^ 2
              )))))
    Xbeta <- X %*% BETA
    xqr <- qr(X)
    Qx <- qr.Q(xqr)
    G <- crossprod(X)
    Ginv <- solve(G)
    svec <- get_s_ldet(G)
    Cmat <- get_Cmat_eigen(X, svec, G, Ginv)
    Y <- Xbeta + rnorm(N)
    XYcp <- crossprod(X, Y)
    abs_XYcp <- abs(XYcp)
    sel <- select_manyXtilde(svec, X, Qx, Y, abs_XYcp, Cmat, Ginv, FDR, offset = offset, nxtil = nxtil)
    c(fdp(sel), fpr(sel), 
      ppv(sel), tpr(sel), 
      length(sel), 1 * one_to_p %in% sel)
  }

res_fmt <- 
  as_tibble(simres) %>%
  mutate(sim_ix = row_number())
res_fmt$pop_Sigma_cnum <- kappa(SigmaGen, exact = TRUE)
res_fmt$k <- k
res_fmt$FDR <- FDR
res_fmt$signal <- magnitude
res_fmt$sigmatype <- sigmatype
res_fmt$rho <- rho
res_fmt$offset <- offset
res_fmt$N <- N
res_fmt$nxtil <- nxtil
res_fmt$p <- p
save(res_fmt, BETA, myseed, file = paste("xtil-", nxtil, "-statcordiff-N", N, 
                                         "-p", p, "-rho", rho,
                                         "-off", offset, '-', sigmatype,
                                         '.RData', 
                                         sep='')
)
