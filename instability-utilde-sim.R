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
nsim <- as.numeric(args[3])
myseed <- as.numeric(args[4])
mycores <- as.numeric(args[5])
sigmatype <- args[6]
k <- as.numeric(args[7])
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
mycores <- max(c(mycores, detectCores() - 1))
options(cores=mycores)
options(mc.cores=mycores)
cat("\n--using "); cat(mycores); cat(" cores\n")

knockoff_select <- function(X, svec, xqr, Y, FDR, offset=1, stat = stat.lasso_lambdadiff, random=TRUE){
  Xtilde <- get_knockoffs_qr(X, svec, xqr = xqr, random = random)
  W <- stat(X, Xtilde, Y)
  return(which(W >= knockoff.threshold(W, FDR, offset=offset)))
}

get_1simfun <- function(N, SigmaGen, BETA, FDR, offset = 1, stat=stat.lasso_coefdiff, random=TRUE){
  k <- sum(abs(BETA) > 0)
  p <- ncol(SigmaGen)
  X <- mvrnorm(N, mu=rep(0, p), Sigma = SigmaGen)
  # FIX X, only vary Utilde
  X <- scale(X, center=T, scale=F)
  X <- scale(X, center=F, scale=apply(X, 2, function(xj) return(sqrt(sum(xj^2)))))
  xqr <- qr(X)
  G <- crossprod(X)
  s_all <- get_all_svec(G)
  Y <- X %*% BETA + rnorm(N)
  fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
  ppv <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))
  tpr <- function(selected) sum(BETA[selected] != 0) / k
  ret <- function(i){
    sel_all <- lapply(s_all, function(svec) return(knockoff_select(X, svec, xqr, Y, FDR, offset=offset, 
                                                                   stat = stat, 
                                                                   random = random)))
    return(sapply(sel_all, function(ss) return(c(fdp = fdp(ss),
                                                 ppv = ppv(ss),
                                                 tpr = tpr(ss),
                                                 nsel = length(ss)))))
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

simres <- vector('list', nsim)
BETA <- vector('numeric', p)
magnitude <- 3.5
BETA[sample(1:p, size = k)] <- magnitude * sample(c(1,-1), size=k, replace=T)
FDR <- 0.1
cat("Nominal FDR = "); cat(FDR)
sim1fun <- get_1simfun(N, SigmaGen, BETA, FDR, offset=0, stat=stat.lasso_coefdiff, random=TRUE)
simres <- mclapply(1:nsim, sim1fun)

res_fmt <- 
  bind_rows(lapply(simres, as_tibble, rownames='meas'),
          .id = 'sim') %>%
  gather(sdp,ldet,equi, key='stype',value='val')%>%
  spread(meas, val)
res_fmt$pop_Sigma_cnum <- kappa(SigmaGen, exact = TRUE)
res_fmt$k <- k
res_fmt$FDR <- FDR
res_fmt$signal <- magnitude
res_fmt$sigmatype <- sigmatype

write.table(res_fmt, file = paste("utilde-unstable-N", N, "-p", p, 
                                  '.csv', 
                                  sep=''),
            quote=F, sep=',', row.names=FALSE)
