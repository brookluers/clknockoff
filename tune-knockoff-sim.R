#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(parallel)
library(knockoff)
library(MASS)
source("tune-knockoff.R")
N <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])
myseed <- as.numeric(args[4])
mycores <- args[5]
sigmatype <- args[6]
betatype <- args[7]
k <- as.numeric(args[8])
r2_betafix <- as.numeric(args[9])
if (is.na(args[10])){
  rhotest <- seq(0.1,0.9,length.out = 8)
} else {
  rhotest <- as.numeric(args[10])
}

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

### Some 'global' values
lambda_ridge <- 0.2
FDR <- 0.2
kindices <- sample(1:p, size = k, replace = FALSE)
k_signs <- sample(c(1,-1), size=k, replace=T)
cat("\nnominal FDR = "); cat(FDR)
if (betatype!='flatfixed'){
  cat("\nsetting beta based on fixed R^2, ")
  cat("\ndesired R^2 = "); cat(r2_betafix)
} 


get_Xgenfunc <- function(SigmaGen, N){
  L <- chol(SigmaGen)
  p <- ncol(SigmaGen)
  f <- function() {
    matrix(rnorm(N*p), nrow=N, ncol=p) %*% L
  }
  return(f)
}

getSigma_exch <- function(p, rho){
  SigmaGen <- matrix(0, nrow=p, ncol=p)
  SigmaGen[lower.tri(SigmaGen, diag = F)] <- rho
  SigmaGen <- SigmaGen + t(SigmaGen) + diag(p)
  return(SigmaGen)
}

getSigma_ar1 <- function(p, rho){
  SigmaGen <- matrix(0, nrow=p, ncol=p)
  SigmaGen <- rho^abs(row(SigmaGen) - col(SigmaGen))
  return(SigmaGen)
}

getSigma_1band <- function(p, rho) {
  Sigma <- matrix(0, nrow=p, ncol=p)
  for (j in 1:(p-1)){
    Sigma[j,j+1] <- rho
  }
  Sigma <- Sigma + t(Sigma)
  diag(Sigma) <- 1
  return(Sigma)
}

getBETA_linear <- function(p, k, beta_max, kindices, k_signs){
  coef_nonzero <- seq(0, beta_max, length.out=k)
  b <- vector('numeric', p)
  b[kindices] <- coef_nonzero * k_signs
  return(b)
}

getBETA_exponential <- function(p, k, beta_max, kindices, k_signs){
  coef_nonzero <- (1-(3/k))^(1:k) * beta_max
  b <- vector('numeric', p)
  b[kindices] <- coef_nonzero * k_signs
  return(b)
}

getBETA_flat <- function(p, k, beta_max, kindices, k_signs){
  b <- vector('numeric', p)
  b[kindices] <- beta_max * k_signs
  return(b)
}

getBETA_linearfixed <- function(p, k, beta_avg, kindices, k_signs){
  b <- vector('numeric', p)
  b[kindices] <- seq(beta_avg / 2, 3 * beta_avg / 2, length.out=k) * k_signs
  return(b)
}

if (sigmatype == 'exch'){
  getSigmaFunc <- getSigma_exch
} else if (sigmatype=='ar1'){
  getSigmaFunc <- getSigma_ar1
} else if (sigmatype == '1band'){
  getSigmaFunc <- getSigma_1band
  max_rho <- min(abs(-1/(2*cos((1:p)*pi/(p+1)))))
  cat("\nbanded SigmaGen: maximum correlation = "); cat(max_rho)
  cat(" to prevent singular SigmaGen\n")
  if(length(rhotest)>1){
    rhotest <- seq(0.1, max_rho * 0.95, length.out=length(rhotest))
  } else if (max(rhotest) > max_rho){
    cat("\ndesired rho is too big, setting to "); cat(max_rho)
    cat(" * 0.95")
    cat("\n")
    rhotest <- max_rho * 0.95
  }
} else {
  cat("unknown Sigma structure\nusing exchangeable\n\n")
  getSigmaFunc <- getSigma_exch
}

if (betatype =='linear'){
  BETAfunc <- getBETA_linear
} else if (betatype=='exponential'){
  BETAfunc <- getBETA_exponential
} else if (betatype=='flat'){
  BETAfunc <- getBETA_flat
} else if (betatype=='constant'){
  BETAfunc <- getBETA_flat
} else  if (betatype=='flatfixed') {
  BETAfunc <- getBETA_flat
  cat("\nfixing |beta_j| = "); cat(r2_betafix);
  cat(" for all j\n")
} else if (betatype =='linearfixed'){
  BETAfunc <- getBETA_linearfixed
  cat("\n average |beta_j| = "); cat(r2_betafix)
  cat("\n |beta_j| between ")
  cat(r2_betafix/2); cat(" and ")
  cat(3*r2_betafix/2)
} else {
  cat("unknown structure for BETA specified\nusing constant magnitudes\n\n")
  BETAfunc <- getBETA_flat
}

if (!(betatype %in% c("flatfixed", "linearfixed"))){
   bmseq <- seq(0.1,20,length.out=200)
   BETA_grid <- sapply(bmseq, function(bmax) return(BETAfunc(p, k, bmax, kindices, k_signs)))
}

onesimrun <- function(SigmaGen, Xgenfunc, BETA, BETA_smaller, BETA_neq_ix, N, FDR, statfunclist) {
  statnames <- names(statfunclist)
  p <- ncol(SigmaGen)
  X <- Xgenfunc()
  X <- knockoff:::normc(X, center=F)
  Sigma <- crossprod(X)
  Xsvd <- knockoff:::canonical_svd(X)
  svec_equi <- rep(min(c(2 * min(Xsvd$d^2), 1))  , p)
  if (p > 100){
    svec_sdp <- create.solve_asdp(Sigma)
  } else{
    svec_sdp <- create.solve_sdp(Sigma)
  }
  ldetGfunc <- get_ldetfun(Sigma)
  ldetGgrad <- get_ldetgrad(Sigma)
  # Optimize over 0 <= s_j <= 1
  Gopt <-
    constrOptim(rep(0.001,p),
                f = ldetGfunc,
                grad = ldetGgrad,
                ui = rbind(diag(p), -diag(p)),
                ci = c(rep(0,p), rep(-1, p)),
                control = list(fnscale = -1))
  Y <- rnorm(N) + X %*% BETA
  sveclist <- setNames(
    list(svec_equi, svec_sdp, Gopt$par),
    c('equi', 'sdp', 'Gdet')
  )
  Xtlist <- lapply(sveclist, function(ss) return(get_knockoffs(ss, X, Xsvd)))
  Wlist <- 
    lapply(statfunclist, function(stf) return(
      lapply(Xtlist, function(Xk) return(stf(X, Xk, Y)))
    ))
  thresh <- 
    lapply(Wlist,
           function(statW) return(
             lapply(statW, function(W) return(knockoff.threshold(W, fdr=FDR, offset=1)))
           ))
  prop_badorder <- 
    lapply(Wlist,
           function(statW) return(
             lapply(statW, function(W) return(
               mean(outer(W, W, FUN = function(a,b) return(a>b))[lower.tri(SigmaGen,diag=F)][BETA_neq_ix] & 
                      BETA_smaller[BETA_neq_ix]) ))
           ))
  sel <- 
    mapply(W_stat = Wlist, thresh_stat = thresh,
           function(W_stat, thresh_stat){
             mapply(W = W_stat, Tval = thresh_stat,
                    function(W, Tval) return(which(W >= Tval)),
                    SIMPLIFY = FALSE)
           }, SIMPLIFY = F
    )
  ret <- lapply(sel,
                function(sel_stat) return(
                  list(
                    fdp = sapply(sel_stat, function(ss) return(fdp(ss))),
                    tpr = sapply(sel_stat, function(ss) return(tpr(ss))),
                    ppv = sapply(sel_stat, function(ss) return(ppv(ss))),
                    nsel = sapply(sel_stat, length),
                    lGdet = sapply(sveclist, function(ss) return(ldetGfunc(ss))),
                    s = do.call('cbind', sveclist)
                  )
                ))
  for(j in seq_along(ret)){
    ret[[j]]$statname <- statnames[j]
    ret[[j]]$prop_badorder <- unlist(prop_badorder[[j]])
  }
  return(ret)
}

SigmaGenList <- lapply(rhotest, function(rho) return(getSigmaFunc(p, rho)))

statfunclist <- setNames(list(stat.olsdiff, 
                              stat.lasso_lambdadiff, 
                              function(X,X_k,y) return(stat.ridge(X,X_k,y,lambda=lambda_ridge))),
                         c('ols','lasso_lambdadiff', 'stat.ridge'))
simparm <- list(N=N, p=p, nsim=nsim, k=k, kindices, k_signs, r2_betafix,
                betatype, sigmatype, myseed=myseed,lambda_ridge,
                statfunclist, mycores=mycores, FDR=FDR)


### Run the simulation for each generative Sigma matrix
for (rj in seq_along(SigmaGenList)){
  SigmaGen <- SigmaGenList[[rj]]
  
  cat("\nupper 5x5 block of SigmaGen: \n")
  print(SigmaGen[1:5,1:5])
  cat("\n")
  
  if (!(betatype %in% c("flatfixed","linearfixed"))) {
     r2_betamax <- apply(BETA_grid, 2, function(x) return(1 - ( 1 / (t(x) %*% SigmaGen %*% x + 1))))
     BETA <- BETA_grid[,which.min(abs(r2_betamax - r2))]
  } else {
     BETA <- BETAfunc(p, k, r2_betafix, kindices, k_signs)
  }
    
  cat("\nBETA = "); cat(paste(round(BETA,3),collapse=', '))
  cat("\nR^2 = "); cat(1 - (1 / (1 + t(BETA) %*% SigmaGen %*% BETA)))
  cat("\n\n")
  Xgenfunc <- get_Xgenfunc(SigmaGen, N)
  
  fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
  ppv <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))
  tpr <- function(selected) sum(BETA[selected] != 0) / k
  
  BETA_smaller <- outer(BETA,BETA,FUN=function(a,b) return(abs(a) < abs(b)))
  BETA_smaller <- BETA_smaller[lower.tri(BETA_smaller,diag=F)]
  eqtol <- 1e-7
  BETA_neq_ix <- outer(BETA, BETA, FUN=function(a,b) return(abs(a -b) > eqtol))
  BETA_neq_ix <- which(BETA_neq_ix[lower.tri(BETA_neq_ix,diag=F)])
  
  res <-
    mclapply(1:nsim, function(i)
      return(
        onesimrun(SigmaGen, Xgenfunc, BETA, BETA_smaller, BETA_neq_ix,
                  N, FDR, statfunclist)
      ))
  save(
    simparm,BETA,
    res,
    SigmaGen,
    file = paste("sim-", sigmatype, "Sigma", "-", betatype, "BETA", "-", 
                 "-rho",rhotest[rj],"-N",N,"-p",p,
                 "-nsim",nsim,".RData",sep = '')
  )
}