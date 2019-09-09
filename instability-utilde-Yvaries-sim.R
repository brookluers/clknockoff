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
k <- sum(abs(BETA) > 0)
p <- ncol(SigmaGen)
X <- mvrnorm(N, mu=rep(0, p), Sigma = SigmaGen)
# FIX X, only vary Utilde
X <- scale(X, center=T, scale=F)
X <- scale(X, center=F, scale=apply(X, 2, function(xj) return(sqrt(sum(xj^2)))))
Xbeta <- X %*% BETA
xqr <- qr(X)
Qx <- qr.Q(xqr)
G <- crossprod(X)
Ginv <- solve(G)
s_all <- 
  list(
    sdp = knockoff::create.solve_sdp(G),
    ldet = get_s_ldet(G)
  )
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
simres <- vector('list', nsim_Y)
for (i in 1:nsim_Y){
  Y <- Xbeta + rnorm(N)
  XYcp <- crossprod(X, Y)
  abs_XYcp <- abs(XYcp)
  simres[[i]] <- 
    foreach(j=1:nsim_givenY, .combine = rbind)  %dopar% {
    res_byS <- matrix(nrow=nstypes, ncol= result_length)
    rownames(res_byS) <- names(s_all)
    colnames(res_byS) <- c('fdp','fpr','ppv','tpr','nsel', selnames, 'U_ix')
    for (j2 in 1:nstypes){
      svec <- s_all[[j2]]
      Cmat <- Cmats_all[[j2]]
      Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
      Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
      Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde
      Utilde <- qr.Q(qr(Utilde))
      Xtilde <- X - X %*% Ginv_S + Utilde %*% Cmat
      W <- as.numeric(abs_XYcp - abs(crossprod(Xtilde, Y))) # simple correlation differences
      sel <- which(W >= knockoff.threshold(W, FDR, offset=offset))
      res_byS[j2,] <- 
        c(fdp(sel), 
          fpr(sel), ppv(sel),
          tpr(sel), length(sel),
          1 * one_to_p %in% sel, j)
    }
    res_byS
  }
}

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
