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
resnames <- c('fdp','fpr','ppv','tpr','nsel', selnames, 'nU')
result_length <- length(resnames)
simres <- vector('list', nsim)
utest <- c(1, 10)
n_utest <- length(utest)
pzeroes <- rep(0, p)
select_avgW <- function(svec, X, Qx, Y, abs_XYcp, Cmat, Ginv, FDR, offset, nU = 1){
  Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
  X_minus_GinvS <- X - X %*% Ginv_S
  Wlist <- vector('list', nU)
  for(j in seq_along(Wlist)){
    Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
    Utilde <- qr.Q(qr(Utilde)) # orthogonalize only, do not project away from proj_xperp(Y)
    Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde, project away from X
    Xtilde <-  X_minus_GinvS + Utilde %*% Cmat
    Wlist[[j]] <- as.numeric(abs_XYcp - abs(crossprod(Xtilde, Y)))
  }
  avgW <- apply(do.call('rbind', Wlist), 2, mean)
  return(which(avgW >= knockoff.threshold(avgW, FDR, offset=offset)))
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
  s_all <-
    list(sdp = knockoff::create.solve_sdp(G),
         ldet = get_s_ldet(G))
  nstypes <- length(s_all)
  Cmats_all <- lapply(s_all, function(svec) {
    return(get_Cmat_eigen(X, svec, G, Ginv))
  })
  Y <- Xbeta + rnorm(N)
  XYcp <- crossprod(X, Y)
  abs_XYcp <- abs(XYcp)
  res_oneY <- matrix(nrow = nstypes * n_utest, ncol = result_length)
  rownames(res_oneY) <- rep(names(s_all), each = n_utest)
  colnames(res_oneY) <- resnames
  for (j in 1:nstypes) {
    svec <- s_all[[j]]
    Cmat <- Cmats_all[[j]]
    sel_by_nU <- lapply(utest, function(uk) return(select_avgW(svec, X, Qx, Y, abs_XYcp, Cmat, Ginv, FDR, offset = offset, nU = uk)))
    res_oneY[((j - 1) * n_utest + 1):(j * n_utest), ] <-
      do.call('rbind', lapply(1:n_utest, function(j2)  
        return(c(fdp(sel_by_nU[[j2]]), fpr(sel_by_nU[[j2]]), ppv(sel_by_nU[[j2]]),
                 tpr(sel_by_nU[[j2]]), length(sel_by_nU[[j2]]), 1 * one_to_p %in% sel_by_nU[[j2]], utest[j2]))))
  }
  res_oneY
}

res_fmt <- 
  as_tibble(simres, rownames='stype')
res_fmt$pop_Sigma_cnum <- kappa(SigmaGen, exact = TRUE)
res_fmt$k <- k
res_fmt$FDR <- FDR
res_fmt$signal <- magnitude
res_fmt$sigmatype <- sigmatype
res_fmt$rho <- rho
res_fmt$offset <- offset
res_fmt$N <- N
res_fmt$p <- p

save(res_fmt, BETA, myseed, file = paste("avgW-statcordiff-N", N, "-p", p, "-rho", rho,
                                         "-off", offset, '-', sigmatype,
                                         '.RData', 
                                         sep='')
)
