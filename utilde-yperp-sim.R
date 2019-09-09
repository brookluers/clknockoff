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
nsim_X <- as.numeric(args[3])
nsim_Y <- as.numeric(args[4])
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
resnames <- c('fdp','fpr','ppv','tpr','nsel', selnames, 'utype', 'y_ix')
result_length <- length(resnames)
utypes <- c('yperp', 'xperp')  #, 'yrpar')
n_utypes <- length(utypes)
simres_byX <- vector('list', nsim_X)
pzeroes <- rep(0, p)
for (i1 in seq_along(simres_byX)){
  X <- mvrnorm(N, mu = pzeroes, Sigma = SigmaGen)
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
  nstypes <- length(s_all)
  Cmats_all <- lapply(s_all, function(svec){
    return(get_Cmat_eigen(X, svec, G, Ginv))
  })
  simres_byX[[i1]] <- foreach(i2 = 1:nsim_Y, .combine=rbind) %dopar% {
    Y <- Xbeta + rnorm(N)
    Yresid <- Y - Qx %*% crossprod(Qx, Y)
    norm_Yresid <- norm(Yresid, type='F')
    Yresid_normed <- Yresid / norm_Yresid
    XYcp <- crossprod(X, Y)
    abs_XYcp <- abs(XYcp)
    res_oneY <- matrix(nrow = nstypes * n_utypes, ncol = result_length)
    rownames(res_oneY) <- rep(names(s_all), each = n_utypes)
    colnames(res_oneY) <- resnames
    for (j2 in 1:nstypes){
      svec <- s_all[[j2]]
      Cmat <- Cmats_all[[j2]]
      Ginv_S <- sweep(Ginv, 2, svec, FUN=`*`)
      X_minus_XGinvS <- X - X %*% Ginv_S
      Utilde <- matrix(rnorm(N*p), nrow=N, ncol=p)
      Utilde <- Utilde - Qx %*% crossprod(Qx, Utilde) #(I - QQ^t) Utilde, project away from X
      Utilde_xperp <- qr.Q(qr(Utilde)) # orthogonalize only, do not project away from proj_xperp(Y)
      ## Choosing Utilde to be orthogonal to proj_X(Y) and proj_Xperp(Y)
      Utilde_yperp <- Utilde - Yresid_normed %*% crossprod(Yresid_normed, Utilde) # project away from proj_xperp(Y)
      #Utilde_yrpar <- cbind(Utilde[,1:(p-1)], Yresid_normed) # Yresid_normed %*% crossprod(Yresid_normed, Utilde)
      #Utilde_yrpar <- Utilde_yrpar - Qx %*% crossprod(Qx, Utilde_yrpar) # parallel to proj_xperp(Y)
      #Utilde_yrpar <- qr.Q(qr(Utilde_yrpar))
      Utilde_yperp <- qr.Q(qr(Utilde_yperp)) # orthogonalize
      Xtilde_yperp <- X_minus_XGinvS + Utilde_yperp %*% Cmat
      Xtilde_xperp <- X_minus_XGinvS + Utilde_xperp %*% Cmat
      #Xtilde_yrpar <- X_minus_XGinvS + Utilde_yrpar %*% Cmat 
      W_yperp <- as.numeric(abs_XYcp - abs(crossprod(Xtilde_yperp, Y))) # simple correlation differences
      W_xperp <- as.numeric(abs_XYcp - abs(crossprod(Xtilde_xperp, Y)))
      #W_yrpar <- as.numeric(abs_XYcp - abs(crossprod(Xtilde_yrpar, Y)))
      sel_yperp <- which(W_yperp >= knockoff.threshold(W_yperp, FDR, offset=offset))
      sel_xperp <- which(W_xperp >= knockoff.threshold(W_xperp, FDR, offset=offset))
      #sel_yrpar <- which(W_yrpar >= knockoff.threshold(W_yrpar, FDR, offset=offset))
      res_oneY[((j2-1) * n_utypes + 1):(j2 * n_utypes),] <- 
        rbind(
          c(fdp(sel_yperp), 
            fpr(sel_yperp), ppv(sel_yperp),
            tpr(sel_yperp), length(sel_yperp),
            1 * one_to_p %in% sel_yperp, 1, i2),
          c(fdp(sel_xperp), 
            fpr(sel_xperp), ppv(sel_xperp),
            tpr(sel_xperp), length(sel_xperp),
            1 * one_to_p %in% sel_xperp, 2, i2))
          #c(fdp(sel_yrpar), 
          #  fpr(sel_yrpar), ppv(sel_yrpar),
          #  tpr(sel_yrpar), length(sel_yrpar),
          #  1 * one_to_p %in% sel_yrpar, 3, i2))
    }
    res_oneY
  }
}

res_fmt <- 
  bind_rows(lapply(simres_byX, as_tibble, rownames='stype'),
          .id = 'x_ix') %>%
  mutate(utype = utypes[utype]) 

res_fmt$pop_Sigma_cnum <- kappa(SigmaGen, exact = TRUE)
res_fmt$k <- k
res_fmt$FDR <- FDR
res_fmt$signal <- magnitude
res_fmt$sigmatype <- sigmatype
res_fmt$rho <- rho
res_fmt$offset <- offset
res_fmt$N <- N
res_fmt$p <- p

save(res_fmt, BETA, myseed, file = paste("utilde-yperp-statcordiff-N", N, "-p", p, "-rho", rho,
                                         "-off", offset, '-',sigmatype,
                                         '.RData', 
                                         sep='')
)
