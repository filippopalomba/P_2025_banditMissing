banditSimul <- function(rounds, A, mu, p, sig.bar, q.bar=0.2, d=0L, rho=0L,
                        regularize=FALSE, oracle.benchmark=FALSE, UCB.type="UCB",
                        light.output=TRUE, sig2.C=2, sig2.R=1) {

  # rounds:           total number of rounds (T)
  # A:                number of actions (|\mathcal{A}|)
  # p:                vector with probability of missingness for each action
  # sig.bar:          variance proxy for rewards
  # q.bar:            lower bound on "propensity score" (proba of not being missing)
  # d:                dimension of the covariate space
  # rho:              scalar governing the (unconditional) correlation bw missingness and rewards
  # regularize:       useful when allowing for missingness and UCB, it prevents the denominator to be 0
  # oracle.benchmark: run oracle algorithm that uses true underlying distributions instead of bonus
  # UCB.type:         chooses the type of UCB algorithm between standard ("UCB") and doubly-robust ("UCB-dr")
  # light.output:     suppresses nuisance estimation from stored output
  # sig2.C:           noise in missingness
  # sig2.R:           noise in rewards
  
  ## Check inputs
  if (length(p) != A) {
    stop("p has to be a list of length A")
  }

  if (!(UCB.type %in% c("UCB", "UCB-dr"))) UCB.type <- "UCB"

  ## Initiate UCB algorithm
  
  # set regularization parameter to avoid ill-defined denominator
  # (only effectively used if 'UCB' + missingness)
  if (any(p != 1)) regularize <- TRUE
  lambda <- 0
  if (isTRUE(regularize)) lambda <- log(rounds)   

  # construct parametrization that guarantees desired rho
  if (rho != 0L) {
    if (any(p == 1)) {
      rho <- 0L
    } else {
      if (d == 0L) d <- 1L
      beta <- unlist(lapply(c(1:A), function(i)
        optimize(calibCorr, interval=c(0, 10), rho.CR=rho, d=d, p=p[i],
                 sig2.C=sig2.C, sig2.R=sig2.R)$minimum))
    }
  }

  ## If UCB-dr is chosen, need to train nuisance estimators on another sample
  if (UCB.type == "UCB-dr") {
    nuisanceEst <- nuisanceExtData(A, mu, sig.bar, p, d, rho, beta, sig2.C, sig2.R)
  } else {
    nuisanceEst <- NULL
  }

  # pull each arm once to initialize algorithm
  params <- prepareAlgParams(rounds=rounds, A=A, mu=mu, sig.bar=sig.bar, q.bar=q.bar, p=p,
                             lambda=lambda, rho=rho, beta=beta, d=d, regularize=regularize,
                             UCB.type=UCB.type, nuisanceEst=nuisanceEst, sig2.C=sig2.C, sig2.R=sig2.R)

  R.t <- params$R.t
  C.t <- params$C.t
  R.til <- params$R.til
  R.hat <- params$R.hat
  P.t <- params$P.t
  N.t <- params$N.t
  lambda <- params$lambda
  bonus <- params$bonus
  mu.hat <- params$mu.hat
  q.hat <- params$q.hat
  beta <- params$beta

  if (isTRUE(oracle.benchmark)) {
    bonus <- oracleBonusGet(rho=rho, N.t=N.t, P.t=P.t, lambda=lambda,
                            beta=beta, A=A, p=p, mu=mu, rounds=rounds, UCB.type=UCB.type)
    R.til <- R.hat + bonus
  }

  R.hat.store <- as.matrix(R.hat)
  R.til.store <- as.matrix(R.til)
  C.store <- C.t
  R.store <- R.t
  R.t.store <- as.list(R.t)
  C.t.store <- as.list(C.t)
  mu.hat.store <- as.list(mu.hat)
  q.hat.store <- as.list(q.hat)
  I.t.store <- c(1:A)
  
  for (t in seq(from=A+1, to=rounds+A, by=1)) {
    
    # play action and observe environment
    I.t <- which.max(R.til)                                         # index of best action in terms of R.tilde
    env <- environmentGet(I=I.t, A=A, mu=mu, sig.bar=sig.bar,
                          p=p, d=d, rho=rho, beta=beta, 
                          sig2.C=sig2.C, sig2.R=sig2.R)             # sample realizations from environment
    
    R.t.store[[I.t]] <- append(R.t.store[[I.t]], env[["R.t"]])      # store rewards
    C.t.store[[I.t]] <- append(C.t.store[[I.t]], env[["C.t"]])      # store missingness mechanism (for each action)
    R.store <- c(R.store, env[["R.t"]])                        # store whole stream of rewards
    C.store <- c(C.store, env[["C.t"]])                        # store whole stream, of missingness mechanism

    if (UCB.type == "UCB-dr") {
      est <- nuisancePredict(env[["X.t"]], nuisanceEst[[I.t]], d)
      q.hat.store[[I.t]] <- append(q.hat.store[[I.t]], est[["q.hat"]])
      mu.hat.store[[I.t]] <- append(mu.hat.store[[I.t]], est[["mu.hat"]])
    }

    # update algorithm parameters for next round
    paramsUpdated <- updateAlgParams(t=t, I.t=I.t, R.hat=R.hat, R.til=R.til, P.t=P.t, N.t=N.t,
                                     R.t.hist=R.t.store[[I.t]], C.t.hist=C.t.store[[I.t]],
                                     mu.hat.hist=mu.hat.store[[I.t]], q.hat.hist=q.hat.store[[I.t]],
                                     rounds=rounds, A=A, lambda=lambda, sig.bar=sig.bar,
                                     q.bar=q.bar, d=d, UCB.type=UCB.type) 
    P.t <- paramsUpdated$P.t
    N.t <- paramsUpdated$N.t
    R.hat <- paramsUpdated$R.hat
    R.til <- paramsUpdated$R.til

    if (isTRUE(oracle.benchmark)) {
      bonus <- oracleBonusGet(rho=rho, N.t=N.t, P.t=P.t, lambda=lambda,
                              beta=beta, A=A, p=p, mu=mu, rounds=rounds, UCB.type=UCB.type)
      R.til <- R.hat + bonus
    }

    # update history
    R.hat.store <- cbind(R.hat.store, R.hat)
    R.til.store <- cbind(R.til.store, R.til)
    I.t.store <- c(I.t.store, I.t)
  }
  
  
  I.t.store <- I.t.store[(A+1):(rounds+A)]
  ## dataset with behavior of each action
  df.store <- data.frame("round" = c(1:(rounds)),
                         "actionChosen" = I.t.store,
                         "R.hat" = t(R.hat.store[,c(2:(rounds+1))]),
                         "R.til" = t(R.til.store[,c(2:(rounds+1))]))
  
  df.long <- df.store %>%
      tidyr::pivot_longer(cols = starts_with("R.hat") | starts_with("R.til"),
                          names_to = c(".value", "action"),
                          names_pattern = "(R\\.hat\\.|R\\.til\\.)(\\d+)")
  
  names(df.long) <- sub("\\.$", "", names(df.long))
  df.long$actionIsChosen <- df.long$actionChosen == df.long$action
  
  df.long$actionIsChosen <- factor(df.long$actionIsChosen, 
                                   levels = c(TRUE, FALSE), 
                                   labels = c("Yes", "No"))

  ## dataset with "observed" behavior
  df.regret <- data.frame("round" = c(1:rounds),
                          "actionChosen" = I.t.store,
                          "censoring" = C.store[(A+1):(rounds+A)],
                          "reward" = R.store[(A+1):(rounds+A)],
                          "regret" = max(mu) - mu[I.t.store]
                          )

  df.regret$regretCumul <- cumsum(df.regret$regret)

  if (rho > 0) {
    rho.sim <- unlist(lapply(c(1:A), function(a) getEmpCorr(a, subset(df.regret, actionChosen == a))))
  } else {
    rho.sim <- 0L
  }
  
  if (isTRUE(light.output)) {
    to.return <- list(df.long=df.long, df.regret=df.regret, actionChosen=I.t.store, rho.sim=rho.sim)
  } else {
    to.return <- list(df.long=df.long, df.regret=df.regret, actionChosen=I.t.store, rho.sim=rho.sim,
                      nuisanceEst=nuisanceEst)
  }
  return(to.return)
}

## function to update UCB parameters at the end of each round
updateAlgParams <- function(t, I.t, P.t, N.t, R.hat, R.til, R.t.hist, C.t.hist, mu.hat.hist, q.hat.hist,
                            rounds, A, lambda, sig.bar, q.bar, d, UCB.type) {
  delta <- 1/rounds^2
  P.t[I.t] <- P.t[I.t] + 1
  N.t[I.t] <- N.t[I.t] + C.t.hist[[length(C.t.hist)]] 
  R.hat[I.t] <- RhatGet(N.t=N.t[I.t], P.t=P.t[I.t], R.t.hist=R.t.hist, C.t.hist=C.t.hist, lambda=lambda,
                        UCB.type=UCB.type, q.hat.hist=q.hat.hist, mu.hat.hist=mu.hat.hist)
  R.til[I.t] <- R.hat[I.t] + bonusGet(rounds=rounds, delta=delta, P.t=P.t[I.t], N.t=N.t[I.t], lambda=lambda,
                                      sig.bar=sig.bar, UCB.type=UCB.type)
  
  return(list(P.t=P.t, N.t=N.t, R.hat=R.hat, R.til=R.til))
}

## Function that samples the random variables of the bandit
environmentGet <- function(I, A, mu, sig.bar, p, d=0, rho=0, beta=NULL, mu.X=0, sig2.C=2, sig2.R=1) {
  # I:        index of the arm being pulled
  # sig.bar:  upper bound on rewards (lower bound is 0)
  # p:        probability of missingness
  # d:        dimension of the covariate space
  # rho:      scalar governing the correlation bw missingness and rewards
  
  ########################################################
  ## Reward and missingness are independent
  if (rho == 0) { 
    R.t <- stats::rnorm(1, mean = mu[I], sd = sig.bar)
    C.t <- stats::rbinom(1, 1, p[I])
    X.t <- NULL
  }
  
  ########################################################
  ## Reward and missingness are dependent
  if (rho > 0) {

    beta <- c(beta[I], rep(0.1, d-1))
    U.C <- rnorm(1, 0, sqrt(sig2.C))
    U.R <- rnorm(1, 0, sqrt(sig2.R))

    X.t <- MASS::mvrnorm(n=1, mu=rep(mu.X, d), Sigma=diag(d))
    W <- sum(X.t * beta)
    sig2.beta <- sum(beta^2)
    
    q.p <- qnorm(1-p[I], mean=0, sd = sqrt(sig2.beta + sig2.C))
    
    C.t <- 1*(W - U.C > q.p)
    R.t <- mu[I] + W + U.R

  }
  
  return(list(R.t=R.t, C.t=C.t, X.t=X.t))
  
}

## function that computes the UCB optimistic bonus
bonusGet <- function(rounds, delta, P.t, N.t, lambda, sig.bar, UCB.type="UCB",
                     K.bar = 1, q.bar=0.2) {

  if (UCB.type == "UCB") {
    b <- (sig.bar / q.bar) * sqrt(2 * log(2/delta) / (N.t + lambda)) + lambda * K.bar / (N.t + lambda)    
  }

  if (UCB.type == "UCB-dr") {
    b <-  (sig.bar / q.bar + sig.bar) * sqrt(2 * log(2/delta) / P.t)  
  }
  
  return(b)
}

## function that initiates the UCB algorithm (just run once)
prepareAlgParams <- function(rounds, A, mu, sig.bar, q.bar, lambda, rho, beta, p=p,
                             d, regularize, UCB.type, nuisanceEst, sig2.C, sig2.R) {

  delta <- 1/rounds^2
  burnin <- lapply(c(1:A), function(ell) environmentGet(I=ell, A=A, mu=mu, sig.bar=sig.bar,
                                                        p=p, d=d, rho=rho, beta=beta, 
                                                        sig2.C=sig2.C, sig2.R=sig2.R))

  ## set regularization parameter and initialize relevant variables
  if (UCB.type == "UCB") {
    
    P.t <- rep(1, A)                             # number of times arm has been pulled (vector)
    C.t <- unlist(lapply(burnin, "[[", "C.t"))   # number of times feedback has been observed (vector)
    N.t <- C.t
    R.t <- unlist(lapply(burnin, "[[", "R.t"))   # mean reward estimate (vector)
    R.hat <- R.t * N.t/ (N.t + lambda)           # apply missingness
    bonus <- bonusGet(rounds=rounds, delta=delta, P.t=P.t, N.t=N.t, lambda=lambda,
                      sig.bar=sig.bar, UCB.type=UCB.type)
    R.til <- R.hat + bonus                       # optimistic mean reward estimate (vector)
    q.hat <- mu.hat <- X.t <- NULL

  } else {
    
    X.t <- lapply(burnin, "[[", "X.t")
    C.t <- unlist(lapply(burnin, "[[", "C.t"))    # number of times feedback has been observed (vector)
    N.t <- C.t
    est <- lapply(c(1:A), function(i) nuisancePredict(X.t[[i]], nuisanceEst[[i]], d))
    mu.hat <- unlist(lapply(est, "[[", "mu.hat")) # get estimated long conditional expectation
    q.hat <- unlist(lapply(est, "[[", "q.hat"))   # get estimated conditional probability of missingness  
    P.t <- rep(1, A)                              # number of times arm has been pulled (vector)
    R.t <- unlist(lapply(burnin, "[[", "R.t"))    # mean reward estimate (vector)
    R.hat <- (C.t * (R.t - mu.hat) / q.hat + mu.hat) / P.t
    bonus <- bonusGet(rounds=rounds, delta=delta, P.t=P.t, N.t=N.t, lambda=lambda,
                      sig.bar=sig.bar, UCB.type=UCB.type)
    R.til <- R.hat + bonus   

  }

  params <- list(R.til = R.til, R.hat = R.hat, P.t = P.t, N.t = N.t, R.t = R.t, C.t = C.t, X.t = X.t,
                 q.hat = q.hat, mu.hat = mu.hat, lambda = lambda, bonus = bonus,
                 beta = beta)

  return(params)
}

RhatGet <- function(N.t, P.t, R.t.hist, C.t.hist, lambda, UCB.type, q.hat.hist=NULL, mu.hat.hist=NULL) {
  if (UCB.type == "UCB") {
    Rhat <- sum(R.t.hist * C.t.hist) / (N.t + lambda)
  }
  
  if (UCB.type == "UCB-dr") {
    Rhat <- sum(C.t.hist * (R.t.hist - mu.hat.hist) / q.hat.hist + mu.hat.hist) / P.t 
  }
  
  return(Rhat)
}

calibCorr <- function(x, rho.CR, d, p, sig2.C=2, sig2.R=1) {
  
  # this code maintains P(C=1) = p and searches for
  # a specific parametrization that yields Corr(C, R) = rho
  
  beta <- c(x, rep(0.1, d))
  sig2.beta <- sum(beta^2)
  q <- qnorm(1-p, mean=0, sd = sqrt(sig2.beta + sig2.C))
  qtil <- q / sqrt(sig2.beta + sig2.C)
  
  rho <- sqrt(sig2.beta / (sig2.beta + sig2.C))
  eps <- rho.CR -  rho * sqrt(sig2.beta / ( (sig2.beta + sig2.R) * p * (1-p))) * dnorm(qtil)
  
  return(eps^2)
}

getEmpCorr <- function(a, df) {
  if (nrow(df) < 100) {
    est <- NA
  } else {
    est <- cor(df$reward, df$censoring)
  }
  return(est)
}

oracleBonusGet <- function(rho, N.t, P.t, lambda, beta, A, p, mu, rounds, UCB.type,
                           sig2.UC = 2, sig2.UR = 1) {
  
  alph <- 1-1/rounds^2
  
  if (rho == 0) {
    
    bonus <- qnorm(alph) * sqrt(sig2.UR) / sqrt(N.t + lambda) 

  } else {

    if (UCB.type == "UCB") {
      sig2.beta <- unlist(lapply(beta, function(b) sum(b^2)))
      sig2.R <- sig2.beta + sig2.UR
      sig2.V <- sig2.beta + sig2.UC
      rho <- sig2.beta / sqrt(sig2.R * sig2.V)
      sig <- (1 - rho^2) * sig2.R
      q.p <- unlist(lapply(c(1:A), function(i) 
        qnorm(1-p[i], mean=0, sd = sqrt(sig2.beta[i] + sig2.UC[i]))))
      bonus <- unlist(lapply(c(1:A), function(i) 
        truncnorm::qtruncnorm(alph, a = q.p[i], mean = 0, sd = sqrt(sig[i])) / sqrt(P.t[i])))

    } else {
      sig2.beta <- unlist(lapply(beta, function(b) sum(b^2)))
      bonus <- qnorm(alph) * sqrt(sig2.beta + sig2.UR) / sqrt(P.t)

    }
  }
  
  return(bonus)
}

nuisanceExtData <- function(A, mu, sig.bar, p, d, rho, beta, sig2.C, sig2.R, aux.data.size=1000) {

  # generate external sample
  aux.data <- list()
  for (a in seq_len(A)) {
    aux <- lapply(c(1:aux.data.size), function(i) environmentGet(I=a, A=A,mu= mu, sig.bar=sig.bar,
                                                                 p=p, d=d, rho=rho, beta=beta, 
                                                                 sig2.C=sig2.C, sig2.R=sig2.R))
    aux <- lapply(aux, function(samp) c(samp[["R.t"]], samp[["C.t"]], samp[["X.t"]]))
    aux <- Reduce(rbind, aux)
    aux.df <- data.frame(aux)
    names(aux.df) <- c("R", "C", paste0("X",c(1:d)))
    aux.data[[a]] <- aux.df
  }

  nuisanceEst <- lapply(aux.data, function(df) nuisanceEstimate(df))  

  return(nuisanceEst)
}

nuisanceEstimate <- function(df) {
  
  # identify all columns starting with "X"
  x_cols <- grep("^X", names(df), value = TRUE)
  formula.mu <- as.formula(paste("R ~", paste(x_cols, collapse = " + ")))
  formula.q <- as.formula(paste("C ~", paste(x_cols, collapse = " + ")))

  # estimate long conditional expectation of rewards using LS
  mu.hat <- lm(formula.mu, data=df)
  
  # estimate conditional probability of missingness using probit
  q.hat <- glm(formula.q, family = binomial(link = "probit"), df)
  
  return(list("mu"=mu.hat, "q"=q.hat))
}

nuisancePredict <- function(X, nuisanceEst, d) {

  X <- as.list(X); names(X) <- paste0("X", c(1:d))
  mu.hat <- predict(nuisanceEst$mu, newdata = X)
  q.hat <- predict(nuisanceEst$q, newdata = X, type = "response")

  return(list(mu.hat=mu.hat, q.hat=q.hat))
}


#### aux funs

muCensGet <- function(rho, d, p, mu, mu.X=0, sig2.C=2, sig2.R=1) {

  beta <- optimize(calibCorr, interval=c(0, 10), rho.CR=rho, d=d, p=p, sig2.C=sig2.C, sig2.R=sig2.C)$minimum
  
  sig2.beta <- sum(beta^2)
  q <- qnorm(1-p, mean=0, sd = sqrt(sig2.beta + sig2.C))
  qtil <- q / sqrt(sig2.beta + sig2.C)
  
  kappa <- sig2.beta / (sig2.beta + sig2.C)
  mu.cens <- mu + kappa * (sqrt(sig2.beta + sig2.C) * dnorm(qtil) / (1-pnorm(qtil)))
  
  return(list(mu.cens=mu.cens, beta=beta))    
}

###################################################################################
# auxiliary function that simulates a draw from a bandit and estimates mean rewards 
# using various estimators
###################################################################################

simulMissing <- function(A, d, p, rho, mu, sims, mu.X=0, sig2.C=2, sig2.R=1,
                           light.output = TRUE) {
  
  mu.cens <- muCensGet(rho, d, p, mu, mu.X, sig2.C, sig2.R)
  beta <- mu.cens$beta
  mu.cens <- mu.cens$mu.cens
  
  draws <- lapply(c(1:sims), function(t) unlist(environmentGet(I=1, A=A, mu=mu, sig.bar=sig.bar,
                                                               p=p, d=d, rho=rho, beta=beta,
                                                               mu.X=mu.X, sig2.C=sig2.C, sig2.R=sig2.R)))
  
  draws <- Reduce(rbind, draws)
  draws <- data.frame(draws)
  draws$R.obs <- draws$R.t * draws$C.t
  names(draws) <- gsub(".t", "", names(draws))
  if (d==1) names(draws)[3] <- "X1"

  nuisanceEst <- nuisanceExtData(A, mu, sig.bar, p, d, rho, beta, sig2.C, sig2.R)[[1]]

  mu.hat <- predict(nuisanceEst$mu, newdata = draws)
  q.hat <- predict(nuisanceEst$q, newdata = draws, type = "response")
  draws$Rtil <- draws$C * (draws$R - mu.hat) / q.hat + mu.hat

  est.oracle <- mean(draws$R)
  est.censor <- mean(subset(draws, C == 1)$R)  
  est.dr <- mean(draws$Rtil)
  
  est.oracle.rolling <- dplyr::cummean(draws$R)
  est.censor.rolling <- cumsum(draws$R*draws$C) / (cumsum(draws$C) + log(sims)) 
  est.dr.rolling <- dplyr::cummean(draws$Rtil)
  
  if (isTRUE(light.output)) {
    return(list(point = c(est.oracle, est.censor, est.dr),
           rolling = cbind(est.oracle.rolling, est.censor.rolling, est.dr.rolling)))
  } else {
    return(list(estimates=c(est.oracle, est.censor, est.dr), df=draws))
  }
}


###########################################################################################
# function that creates and saves .tex and .png files of the figures illustrating
# the cumulative regret and probability of choosing the best action
###########################################################################################

graphsCreate <- function(res, path.fig, spec) {
  
  lw.par <- 5
  text.size <- 20
  xlab.text <- "round $(t)$"
  ylab.text.regr <- "regret"
  ylab.text.prob <- "$\\mathbb{P}[A_t = a^\\star]$"
  legend.ucb.label <- "$\\pi=\\pi^{\\mathsf{UCB}}\\:\\:\\:$"
  legend.ucbdr.label <- "$\\pi=\\pi^{\\mathsf{DR}}\\:\\:\\:$"
  legend.oracle.label <- "$\\pi=\\pi_{\\star}^{\\mathsf{UCB}}\\:\\:\\:$"
  legend.super.oracle.label <- "$\\pi=\\pi_{\\star}^{\\mathsf{DR}}\\:\\:\\:$"
  legend.title <- ""
  n.labs <- 2L
  nrows.legend <- 1L
  my.palette <- wes_palette(n=2L, name="GrandBudapest1")
  
  title.m.x <- 10L # x axis title margin from labels
  title.m.y <- 15L # y axis title margin from labels
  text.m.x <- 10L # x axis label margin from ticks
  text.m.y <- 5L # x axis label margin from ticks
  
  leg.sp.x <- 15.0 # horizontal space between legend keys
  
  if (!(spec %in% c("noCen", "indepCen", "depCen"))) stop("spec %in% c('noCen', 'indepCen', 'depCen')")
  
  res.UCB <- res[[paste0(spec, "_UCB")]]
  res.Oracle <- res[[paste0(spec, "_Oracle")]]
  if (spec == "depCen") {
    res.UCBdr <- res[[paste0(spec, "_UCBdr")]]
    res.UCBdr2 <- res[[paste0(spec, "_OracleUCBdr")]]
    my.palette <- wes_palette(n=4L, name="GrandBudapest1")[c(1,3,4)]
  } else {
    res.UCBdr <- NULL
  }
  
  rounds <- nrow(res.UCB[[1]]$df.regret)

  # cumulative regret
  res.regret <- lapply(res.UCB, function(aux) aux$df.regret$regretCumul)
  res.regret <- colMeans(Reduce(rbind, res.regret))
  df.regret <- data.frame("round" = c(1:rounds), "regret" = res.regret, type = "UCB")
  
  res.regret <- lapply(res.Oracle, function(aux) aux$df.regret$regretCumul)
  res.regret <- colMeans(Reduce(rbind, res.regret))
  df.regret2 <- data.frame("round" = c(1:rounds), "regret" = res.regret, type = "Oracle")
  
  toplot <- rbind(df.regret, df.regret2)
  toplot$type <- factor(ifelse(toplot$type == "UCB", 0, 1),
                        levels = c(0, 1),
                        labels = c(legend.ucb.label, legend.oracle.label))  

  if (is.null(res.UCBdr) == FALSE) {
    res.regret <- lapply(res.UCBdr, function(aux) aux$df.regret$regretCumul)
    res.regret <- colMeans(Reduce(rbind, res.regret))
    df.regret3 <- data.frame("round" = c(1:rounds), "regret" = res.regret, type = "UCB-DR")

    res.regret <- lapply(res.UCBdr2, function(aux) aux$df.regret$regretCumul)
    res.regret <- colMeans(Reduce(rbind, res.regret))
    df.regret4 <- data.frame("round" = c(1:rounds), "regret" = res.regret, type = "UCB-DR Oracle")
    
    toplot <- rbind(df.regret, df.regret3, df.regret4)
    toplot$type <- factor(ifelse(toplot$type == "UCB", 0, ifelse(toplot$type == "UCB-DR", 1, 2)),
                          levels = c(0, 1, 2),
                          labels = c(legend.ucb.label, legend.ucbdr.label, legend.super.oracle.label))
    n.labs <- 4L
    nrows.legend <- 2L
  }

  
  p <- ggplot(toplot) +
    geom_line(aes(x=round, y=regret, group=type, color=type), linewidth = lw.par) +
    xlab(xlab.text) + ylab(ylab.text.regr) +
    scale_color_manual(values=my.palette, name=legend.title) +
    guides(colour = guide_legend(nrow=nrows.legend, byrow=TRUE)) +
    theme(legend.position="bottom",
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size, margin = margin(t = title.m.x)),
          axis.title.y = element_text(size = text.size, margin = margin(r = title.m.y)),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(margin = margin(t = text.m.x)),
          axis.text.y = element_text(margin = margin(r = text.m.y)),
          legend.title = element_text(size=text.size),
          legend.text = element_text(size=text.size),
          legend.spacing.x = unit(leg.sp.x, "cm"))
  
  ggsave(paste0(path.fig, "cumulRegret_", spec, ".png"), plot=p)
  
  tikz(file = paste0(path.fig, "cumulRegret_", spec, ".tex"), width = 7, height = 7)
  plot(p)
  dev.off()

  # zoom-in
  if (is.null(res.UCBdr) == FALSE) {

    toplot <- subset(toplot, type %in% c(legend.ucbdr.label, legend.super.oracle.label))
    
    p <- ggplot(toplot) +
      geom_line(aes(x=round, y=regret, group=type, color=type), linewidth = lw.par) +
      xlab(xlab.text) + ylab(ylab.text.regr) +
      scale_color_manual(values=my.palette[2:3], name=legend.title) +
      guides(colour = guide_legend(nrow=1L, byrow=TRUE)) +
      theme(legend.position="bottom",
            axis.text=element_text(size=text.size),
            axis.title.x = element_text(size = text.size, margin = margin(t = title.m.x)),
            axis.title.y = element_text(size = text.size, margin = margin(r = title.m.y)),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(margin = margin(t = text.m.x)),
            axis.text.y = element_text(margin = margin(r = text.m.y)),
            legend.title = element_text(size=text.size),
            legend.text = element_text(size=text.size),
            legend.spacing.x = unit(leg.sp.x, "cm"))
    
    ggsave(paste0(path.fig, "cumulRegret_", spec, "_zoom.png"), plot=p)
    
    tikz(file = paste0(path.fig, "cumulRegret_", spec, "_zoom.tex"), width = 7, height = 7)
    plot(p)
    dev.off()
    
  }
  
  res.action <- lapply(res.UCB, function(aux) aux$actionChosen == which.max(mu))
  res.action <- colMeans(Reduce(rbind, res.action))
  df.action <- data.frame("round" = c(1:rounds), "action" = res.action, type = "UCB")
  
  res.action <- lapply(res.Oracle, function(aux) aux$actionChosen == which.max(mu))
  res.action <- colMeans(Reduce(rbind, res.action))
  df.action2 <- data.frame("round" = c(1:rounds), "action" = res.action, type = "Oracle")
  
  toplot <- rbind(df.action, df.action2)
  toplot$type <- factor(ifelse(toplot$type == "UCB", 0, 1),
                        levels = c(0, 1),
                        labels = c(legend.ucb.label, legend.oracle.label))

  if (is.null(res.UCBdr) == FALSE) {
    res.action <- lapply(res.UCBdr, function(aux) aux$actionChosen == which.max(mu))
    res.action <- colMeans(Reduce(rbind, res.action))
    df.action3 <- data.frame("round" = c(1:rounds), "action" = res.action, type = "UCB-DR")

    res.action <- lapply(res.UCBdr2, function(aux) aux$actionChosen == which.max(mu))
    res.action <- colMeans(Reduce(rbind, res.action))
    df.action4 <- data.frame("round" = c(1:rounds), "action" = res.action, type = "UCB-DR Oracle")
    
    toplot <- rbind(df.action, df.action3, df.action4)
    toplot$type <- factor(ifelse(toplot$type == "UCB", 0, ifelse(toplot$type == "UCB-DR", 1, 2)),
                          levels = c(0, 1, 2),
                          labels = c(legend.ucb.label, legend.ucbdr.label, legend.super.oracle.label))
  }

  p <- ggplot(toplot) +
    geom_line(aes(x=round, y=action, group=type, color=type)) +
    #geom_line(aes(x=round, y=action, group=type, color=type), alpha=0.4) +
    #geom_smooth(aes(x=round, y=action, group=type, color=type), method='lm',
    #            formula= y~log(x), show.legend=FALSE, linewidth = lw.par) +
    geom_hline(aes(yintercept=1), alpha=0.4) +
    xlab(xlab.text) + ylab(ylab.text.prob) +
    scale_color_manual(values=my.palette, name=legend.title) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth=lw.par), nrow=nrows.legend, byrow=TRUE)) +
    theme(legend.position="bottom",
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size, margin = margin(t = title.m.x)),
          axis.title.y = element_text(size = text.size, margin = margin(r = title.m.y)),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(margin = margin(t = text.m.x)),
          axis.text.y = element_text(margin = margin(r = text.m.y)),
          legend.title = element_text(size=text.size),
          legend.text = element_text(size=text.size),
          legend.spacing.x = unit(leg.sp.x, "cm"))
  ggsave(paste0(path.fig, "probBestAction_", spec, ".png"), plot=p)
  
  tikz(file = paste0(path.fig, "probBestAction_", spec, ".tex"), width = 7, height = 7)
  plot(p)
  dev.off()
  
  # zoom-in
  if (is.null(res.UCBdr) == FALSE) {
    toplot <- subset(toplot, type %in% c(legend.ucbdr.label, legend.super.oracle.label))
    p <- ggplot(toplot) +
      geom_line(aes(x=round, y=action, group=type, color=type)) +
      #geom_line(aes(x=round, y=action, group=type, color=type), alpha=0.4) +
      #geom_smooth(aes(x=round, y=action, group=type, color=type), method='lm',
      #            formula= y~log(x), show.legend=FALSE, linewidth = lw.par) +
      geom_hline(aes(yintercept=1), alpha=0.4) +
      xlab(xlab.text) + ylab(ylab.text.prob) +
      scale_color_manual(values=my.palette[2:3], name=legend.title) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth=lw.par), nrow=1L, byrow=TRUE)) +
      theme(legend.position="bottom",
            axis.text=element_text(size=text.size),
            axis.title.x = element_text(size = text.size, margin = margin(t = title.m.x)),
            axis.title.y = element_text(size = text.size, margin = margin(r = title.m.y)),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(margin = margin(t = text.m.x)),
            axis.text.y = element_text(margin = margin(r = text.m.y)),
            legend.title = element_text(size=text.size),
            legend.text = element_text(size=text.size),
            legend.spacing.x = unit(leg.sp.x, "cm"))
    ggsave(paste0(path.fig, "probBestAction_", spec, "_zoom.png"), plot=p)
    
    tikz(file = paste0(path.fig, "probBestAction_", spec, "_zoom.tex"), width = 7, height = 7)
    plot(p)
    dev.off()
  }
}

###########################################################################################
# function that creates and saves .tex and .png files of the figures illustrating
# the 
###########################################################################################


intGet <- function(mat) {
  mu <- as.matrix(rowMeans(mat))
  bounds <- t(apply(mat, 1, quantile, probs = c(0.025, 0.975)))
  cbind(mu, bounds)
}


graphsCreate2 <- function(aux.1, aux.2, plims) {
  lab.oracle <- "$\\check{R}_a(t)\\:\\:$"
  lab.naive <- "$\\widehat{R}_a(t)\\:\\:$"
  lab.dr <- "$\\widehat{R}_a^{\\mathsf{DR}}(t)\\:\\:$"
  
  df <- data.frame(round = c(1:iters),
                   estimator = c(rep(lab.oracle, iters), rep(lab.naive, iters), rep(lab.dr, iters)),
                   estimate = aux.1,
                   action = "$a=1$")
  
  df.2 <- data.frame(round = c(1:iters),
                     estimator = c(rep(lab.oracle, iters), rep(lab.naive, iters), rep(lab.dr, iters)),
                     estimate = aux.2,
                     action = "$a=2$")
  df <- rbind(df, df.2)
  
  colnames(df) <- c("round", "estimator", "mean", "lb", "ub", "action")
  
  df.sub <- subset(df, round > 20)
  
  lw.par <- 5
  text.size <- 13
  title.m.x <- 10L # x axis title margin from labels
  title.m.y <- 15L # y axis title margin from labels
  text.m.x <- 10L # x axis label margin from ticks
  text.m.y <- 5L # x axis label margin from ticks
  
  leg.sp.x <- 15.0 # horizontal space between legend keys
  
  my.palette <- c(wes_palette(n=3L, name="BottleRocket2")[2],
                  wes_palette(n=3L, name="BottleRocket2")[3],
                  wes_palette(n=3L, name="BottleRocket2")[1])
  
  mu.df <- data.frame(action = c("$a=1$", "$a=2$"),
                      intercept = c(0.5, 1))
  plim.df <- data.frame(action = c("$a=1$", "$a=2$"),
                        intercept = plims)
  
  p <- ggplot(df.sub, aes(x = round, y = mean, color = estimator, fill = estimator, group = estimator)) +
    geom_hline(
      data    = mu.df,
      aes(yintercept = intercept, linetype = action),
      linetype="dotted") +
    facet_wrap(~action) +
    geom_ribbon(aes(ymin = lb, ymax = ub),
                alpha = 0.3, colour=NA) +
    geom_line(linetype = "dotted") +
    scale_y_continuous(
      breaks = c(-1, 0, 0.5, 1, 1.5),
      labels = c("-1", "0", "$\\theta_1$", "$\\theta_2$", 1.5)
    ) +
    scale_fill_manual(values=my.palette, name="") +
    scale_color_manual(values=my.palette, name="") +
    xlab("round $(t)$") + ylab("estimate") +
    theme(legend.position="bottom",
          strip.background = element_blank(),
          strip.text = element_text(size = text.size),
          axis.text=element_text(size=text.size),
          axis.title.x = element_text(size = text.size, margin = margin(t = title.m.x)),
          axis.title.y = element_text(size = text.size, margin = margin(r = title.m.y)),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(margin = margin(t = text.m.x), angle = 45),
          axis.text.y = element_text(margin = margin(r = text.m.y)),
          legend.title = element_text(size=text.size),
          legend.text = element_text(size=text.size),
          legend.spacing.x = unit(leg.sp.x, "cm"))
  
  ggsave(paste0(path.fig, "estimatorsWorm.png"), plot=p)
  
  tikz(file = paste0(path.fig, "estimatorsWorm.tex"), width = 7, height = 5)
  plot(p)
  dev.off()
}
