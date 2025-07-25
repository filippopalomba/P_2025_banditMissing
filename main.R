###############################################################################
## Author: Filippo Palomba
## Last modified: 2025-07-25
## Script for simulating MABs with missing rewards
###############################################################################
rm(list=ls(all=TRUE))

##########################################
# Load stuff 
require(pacman)
pacman::p_load(ggplot2, tidyverse, parallel, truncnorm, wesanderson, tikzDevice,
               reshape2)
RNGkind("L'Ecuyer-CMRG")
theme_set(theme_bw())

options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))

##########################################
# Set paths
path <- "YOUR_PATH"
path.code <- paste0(path, "code/")
cores <- 4L
sims <-500L
rounds <- 1000L
analysisDO <- TRUE
plotDO <- TRUE

path.fig  <- paste0(path, "figures/")
path.out  <- paste0(path, "tables/")

source(paste0(path.code, "funs.R"))
set.seed(8894)

##########################################
# Run simulation

A <- 2L                    # number of actions
d <- 1L                    # number of covariates
sig.bar <- 1/2             # variance proxy for sG rewards
mu <- c(0.5, 1)
p.nocen <- rep(1, A)
p.cen <- c(0.2, 0.95)
rho <- 0.2


res <- list()

if (isTRUE(analysisDO)) {
  
  ############################################
  ## No missing data
  ############################################
  cat("No Missing Data - UCB\n")
  res[["noCen_UCB"]] <- parallel::mclapply(c(1:sims), 
                                           function(i) 
                                             banditSimul(rounds=rounds, A=A, mu=mu, p=p.nocen, sig.bar=sig.bar, 
                                                         regularize=FALSE),
                                           mc.cores = cores, mc.set.seed = TRUE
                                           )
  
  aux <- res[["noCen_UCB"]]
  
  cat("No Missing Data - Oracle UCB\n")
  res[["noCen_Oracle"]] <- parallel::mclapply(c(1:sims), 
                                              function(i) 
                                                banditSimul(rounds=rounds, A=A, mu=mu, p=p.nocen, sig.bar=sig.bar, 
                                                            regularize=FALSE, oracle.benchmark = TRUE),
                                              mc.cores = cores, mc.set.seed = TRUE
                                              )
  
  ############################################
  ## Reward-independent missingness
  ############################################

  cat("Reward-Independent Missingness - UCB\n")
  res[["indepCen_UCB"]] <- parallel::mclapply(c(1:sims), 
                                              function(i) 
                                                banditSimul(rounds=rounds, A=A, mu=mu, p=p.cen, 
                                                            sig.bar=sig.bar, regularize=TRUE),
                                              mc.cores = cores, mc.set.seed = TRUE
                                              )
  
  cat("Reward-Independent Missingness - Oracle UCB\n")
  res[["indepCen_Oracle"]] <- parallel::mclapply(c(1:sims), 
                                                 function(i) 
                                                   banditSimul(rounds=rounds, A=A, mu=mu, p=p.cen, sig.bar=sig.bar, 
                                                               regularize=TRUE, oracle.benchmark = TRUE),
                                                 mc.cores = cores, mc.set.seed = TRUE
                                                 )
  
  ############################################
  ## Reward-dependent missingness
  ############################################

  cat("Reward-Dependent Missingness - UCB\n")
  res[["depCen_UCB"]] <- parallel::mclapply(c(1:sims),
                                               function(i) 
                                                 banditSimul(rounds=rounds, A=A, mu=mu, p=p.cen, rho=rho, d=d,
                                                             sig.bar=sig.bar, regularize=TRUE),
                                               mc.cores = cores, mc.set.seed = TRUE
                                            )
  
  cat("Reward-Dependent Missingness - Oracle UCB\n") # not used
  res[["depCen_Oracle"]] <- parallel::mclapply(c(1:sims),
                                               function(i)
                                                 banditSimul(rounds=rounds, A=A, mu=mu, p=p.cen, rho=rho, d=d,
                                                             sig.bar=sig.bar, UCB.type="UCB",
                                                             oracle.benchmark = TRUE),
                                               mc.cores = cores, mc.set.seed = TRUE
                                               )
  
  cat("Reward-Dependent Missingness - DR UCB\n")
  res[["depCen_UCBdr"]] <- parallel::mclapply(c(1:sims),
                                              function(i) 
                                                banditSimul(rounds=rounds, A=A, mu=mu, p=p.cen, rho=rho, d=d,
                                                            sig.bar=sig.bar, UCB.type="UCB-dr"),
                                              mc.cores = cores, mc.set.seed = TRUE
                                              )
  
  cat("Reward-Dependent Missingness - Oracle DR UCB\n")
  res[["depCen_OracleUCBdr"]] <- parallel::mclapply(c(1:sims),
                                              function(i) 
                                                banditSimul(rounds=rounds, A=A, mu=mu, p=p.cen, rho=rho, d=d,
                                                            sig.bar=sig.bar, UCB.type="UCB-dr",
                                                            oracle.benchmark = TRUE),
                                              mc.cores = cores, mc.set.seed = TRUE
                                              )
  
  save(res, file=paste0(path.out,"/simRes.RData"))
}

if (isTRUE(plotDO)) {

  load(paste0(path.out, "simRes.RData"))
  
  graphsCreate(spec="noCen",
               res=res,
               path.fig=path.fig)
  
  graphsCreate(spec="indepCen",
               res=res,
               path.fig=path.fig)

  graphsCreate(spec="depCen",
               res=res,
               path.fig=path.fig)
}


set.seed(8894)

res <- lapply(c(1:sims), function(i) 
  simulMissing(A=1, d=d, p=p.cen[1], rho=rho, mu=mu[1], sims=rounds)$rolling)
res.1 <- Reduce(cbind, lapply(res, function(x) x[, 1]))
res.2 <- Reduce(cbind, lapply(res, function(x) x[, 2]))
res.3 <- Reduce(cbind, lapply(res, function(x) x[, 3]))
aux.1 <- rbind(intGet(res.1), intGet(res.2), intGet(res.3))

res <- lapply(c(1:sims), function(i) 
  simulMissing(A=1, d=d, p=p.cen[2], rho=rho, mu=mu[2], sims=rounds)$rolling)
res.1 <- Reduce(cbind, lapply(res, function(x) x[, 1]))
res.2 <- Reduce(cbind, lapply(res, function(x) x[, 2]))
res.3 <- Reduce(cbind, lapply(res, function(x) x[, 3]))
aux.2 <- rbind(intGet(res.1), intGet(res.2), intGet(res.3))

plims <- c(muCensGet(rho, d, p=p.cen[1], mu=mu[1])$mu.cens,
           muCensGet(rho, d, p=p.cen[2], mu=mu[2])$mu.cens)

graphsCreate2(aux.1, aux.2, plims)

currentfiles <- list.files(path.fig, ".tex$")
setwd(path.fig)
