#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(parallel)
source("opt-G.R")
N <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])
myseed <- as.numeric(args[4])
mycores <- args[5]
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

k <- min(c(floor(p/1.5), 30))
A <- 3.5
FDR <- 0.2
cat("\nFDR = "); cat(FDR);
kcols <- sample(1:p, size = k, replace = FALSE)
BETA <- vector('numeric', length=p)
BETA[kcols] <- sample(c(A, -A), size = k, replace = TRUE)

fdp <- function(selected) sum(BETA[selected] == 0) / max(1, length(selected))
ppv <- function(selected) sum(BETA[selected] != 0) / max(1, length(selected))
tpr <- function(selected) sum(BETA[selected] != 0) / k

rhotest <- seq(0.1,0.9,length.out = 8)
SigmaGenList <- lapply(rhotest, function(rho) {
  SigmaGen <- matrix(0, nrow=p, ncol=p)
  SigmaGen <- rho^abs(row(SigmaGen) - col(SigmaGen))
  return(SigmaGen)
})
# top eigenvalue of Sigma: 1 + (p-1) * a
# other eigenvalues of Sigma: 1 - a 
simparm <- list(N=N, p=p, BETA=BETA, nsim=nsim, k=k, myseed=myseed,
                mycores=mycores, A=A, FDR=FDR)
for (rj in seq_along(SigmaGenList)){
  SigmaGen <- SigmaGenList[[rj]]
  res <- mclapply(1:nsim, function(i) return(onesimrun(SigmaGen, BETA, N, FDR, statfunc=stat.olsdiff,statname='ols')))
  res <- c(res, 
           mclapply(1:nsim, function(i) return(onesimrun(SigmaGen,BETA,N,FDR, statfunc=stat.lasso_lambdadiff,statname='lasso_lambdadiff'))))
  save(simparm, res, SigmaGen,
       file = paste("sim-opt-G-ar1-rho", rhotest[rj],
                    "-N", N, "-p", p, '-nsim', nsim, '.RData', sep=''))
}