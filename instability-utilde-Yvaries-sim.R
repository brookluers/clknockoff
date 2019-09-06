#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(parallel)
library(knockoff)
library(MASS)
source("tune-knockoff.R")
source("pop-cov-funcs.R")
N <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim_Y <- as.numeric(args[3])
nsim_givenY <- as.numeric(args[4])
myseed <- as.numeric(args[5])
mycores <- as.numeric(args[6])
sigmatype <- args[7]
k <- as.numeric(args[8])
offset <- as.numeric(args[9])
if (is.na(args[10])){
  rho <- 0.8
} else {
  rho <- as.numeric(args[10])
}
RNGkind("L'Ecuyer-CMRG")
set.seed(myseed)
cat("\nN = "); cat(N)
cat("\np = "); cat(p)
cat("\nnsim_Y = "); cat(nsim_Y)
cat("\nseed = "); cat(myseed)
cat("\ndesired cores = "); cat(mycores);
cat("\ndetectCores() = "); cat(detectCores())
mycores <- min(c(mycores, detectCores() - 1))
options(cores=mycores)
options(mc.cores=mycores)
cat("\n--using "); cat(mycores); cat(" cores\n")


get_1simfun <- function(N, SigmaGen, BETA, FDR, nsim_givenY = 500, offset = 1, random=TRUE){
  k <- sum(abs(BETA) > 0)
  p <- ncol(SigmaGen)
  X <- mvrnorm(N, mu=rep(0, p), Sigma = SigmaGen)
  # FIX X, only vary Utilde
  X <- scale(X, center=T, scale=F)
  X <- scale(X, center=F, scale=apply(X, 2, function(xj) return(sqrt(sum(xj^2)))))
  Xbeta <- X %*% BETA
  xqr <- qr(X)
  G <- crossprod(X)
  Ginv <- solve(G)
  s_all <- get_all_svec(G)
  Cmats_all <- lapply(s_all, function(svec){
    return(get_Cmat_eigen(X, svec, G, Ginv))
  })
  nstypes <- length(s_all)
  fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
  ppv <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))
  tpr <- function(selected) sum(BETA[selected] != 0) / k
  fpr <- function(selected) sum(BETA[selected] == 0) / (p - k)
  result_length <- 6 + p
  one_to_p <- 1:p
  selnames <- paste('sel', 1:p, sep='')
  ret <- function(i) {
    Y <- Xbeta + rnorm(N)
    XYcp <- crossprod(X, Y)
    abs_XYcp <- abs(XYcp)
    res_given_Y <- matrix(nrow = nsim_givenY * nstypes, ncol = result_length)
    rownames(res_given_Y) <- rep(names(s_all), nsim_givenY)
    colnames(res_given_Y) <- c('fdp','fpr','ppv','tpr','nsel',selnames,'U_ix')
    for (j in 1:nsim_givenY){
      sel_all <- mapply(svec = s_all, Cmat = Cmats_all, 
                        FUN = function(svec, Cmat) {
          Xtilde <- get_knockoffs_qr(X, svec, xqr = xqr, random = random, Cmat = Cmat, Ginv = Ginv, G=G)
          W <- as.numeric(abs_XYcp - abs(crossprod(Xtilde, Y))) # simple correlation differences
          return(which(W >= knockoff.threshold(W, FDR, offset=offset)))
        }, SIMPLIFY = FALSE)
      for (j2 in 1:nstypes){
        res_given_Y[nstypes * (j - 1) + j2, ] <- 
          c(fdp(sel_all[[j2]]), 
            fpr(sel_all[[j2]]), ppv(sel_all[[j2]]),
            tpr(sel_all[[j2]]), length(sel_all[[j2]]),
            1 * one_to_p %in% sel_all[[j2]], j)
      }
    }
    return(res_given_Y)
  }
  return(ret)
}

if (sigmatype == 'exch'){
  SigmaGen <- get_exch(p, rho)
} else if (sigmatype=='ar1'){
  SigmaGen <- get_ar(p, rho)
} else {
  cat("unknown Sigma structure\nusing exchangeable\n\n")
  SigmaGen <- get_exch(p, rho)
}

simres <- vector('list', nsim_Y)
BETA <- vector('numeric', p)
magnitude <- 3.5
BETA[sample(1:p, size = k)] <- magnitude * sample(c(1,-1), size=k, replace=T)
FDR <- 0.1
cat("Nominal FDR = "); cat(FDR); cat("\n")
sim1fun <- get_1simfun(N, SigmaGen, BETA, FDR, nsim_givenY = 50, 
                       offset=offset, random=TRUE)
simres <- mclapply(1:nsim_Y, sim1fun)
#lapply(1:nsim_Y, sim1fun)

res_fmt <- 
  bind_rows(lapply(simres, as_tibble, rownames='stype'),
            .id = 'Y_ix')

res_fmt$pop_Sigma_cnum <- kappa(SigmaGen, exact = TRUE)
res_fmt$k <- k
res_fmt$FDR <- FDR
res_fmt$signal <- magnitude
res_fmt$sigmatype <- sigmatype
res_fmt$rho <- rho
res_fmt$offset <- offset
res_fmt$N <- N
res_fmt$p <- p

save(res_fmt, BETA, myseed, file = paste("utilde-varyY-statcordiff-N", N, "-p", p, "-rho", rho,
                                 "-off", offset, '-',sigmatype,
                                 '.RData', 
                                 sep='')
)
