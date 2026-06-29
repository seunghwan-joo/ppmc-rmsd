mg.twopl.mcmc <- function(data, group,
                          it = 1000,
                          burn = 500,
                          chains = 3,
                          cand_th = .75,
                          cand_mn = .02,
                          cand_std = .02,
                          cand_a = .15,
                          cand_b = .25){
  
  X <- as.matrix(data)
  N <- nrow(X)
  J <- ncol(X)
  
  Grp <- as.integer(factor(group))
  ngrp <- length(unique(Grp))
  
  X_na <- is.na(X)
  
  # ---- Fast 2PL + log-likelihood (handles missing) ----
  twopl_fast <- function(a, b, th) {
    # eta = a*(th - b) = th*a - a*b
    eta <- tcrossprod(th, a) - matrix(a * b, nrow = length(th), ncol = length(a), byrow = TRUE)
    plogis(eta)
  }
  
  loglike_mat_fast <- function(a, b, th) {
    p <- twopl_fast(a, b, th)
    ll <- X * log(p) + (1 - X) * log1p(-p)
    if (any(X_na)) ll[X_na] <- 0
    ll
  }
  
  # ---- Priors (log) ----
  logprior_alpha <- function(alpha) log(dbeta((alpha - .25) / (3 - .25), 1.5, 1.5))
  logprior_beta  <- function(beta)  log(dbeta((beta + 4) / (4 + 4), 2, 2))
  
  initiate <- function() {
    list(
      alpha = rep(1, J),
      beta  = rep(0, J),
      theta = rnorm(N),
      mn    = rep(0, ngrp),
      std   = rep(1, ngrp)
    )
  }
  
  drawtheta_fast <- function(alpha0, beta0, theta0, mn0, std0, acpttheta) {
    theta1 <- theta0 + cand_th * rnorm(N)
    
    thmu <- mn0[Grp]
    thsd <- std0[Grp]
    
    d0 <- rowSums(loglike_mat_fast(alpha0, beta0, theta0)) +
      dnorm(theta0, mean = thmu, sd = thsd, log = TRUE)
    
    d1 <- rowSums(loglike_mat_fast(alpha0, beta0, theta1)) +
      dnorm(theta1, mean = thmu, sd = thsd, log = TRUE)
    
    acc <- d1 - d0
    accind <- (log(runif(N)) < acc)
    
    theta0[accind] <- theta1[accind]
    acpttheta[accind] <- acpttheta[accind] + 1L
    
    list(theta = theta0, acptth = acpttheta)
  }
  
  drawlambda_fast <- function(theta0, mn0, std0, acptmn, acptstd) {
    mn1  <- mn0  + cand_mn  * rnorm(ngrp)
    std1 <- std0 + cand_std * rnorm(ngrp)
    
    # identification
    mn1[1]  <- 0
    std1[1] <- 1
    
    thmu0 <- mn0[Grp];  thsd0 <- std0[Grp]
    thmu1 <- mn1[Grp];  thsd1 <- std1[Grp]
    
    # update mn
    d0 <- sum(dnorm(theta0, mean = thmu0, sd = thsd0, log = TRUE))
    d1 <- sum(dnorm(theta0, mean = thmu1, sd = thsd1, log = TRUE))
    if (log(runif(1)) < (d1 - d0)) {
      mn0 <- mn1
      acptmn <- acptmn + 1L
      thmu0 <- mn0[Grp]  # refresh if accepted
    }
    
    # update std
    d0 <- sum(dnorm(theta0, mean = thmu0, sd = thsd0, log = TRUE))
    d1 <- sum(dnorm(theta0, mean = thmu0, sd = thsd1, log = TRUE))
    if (log(runif(1)) < (d1 - d0)) {
      std0 <- std1
      acptstd <- acptstd + 1L
    }
    
    list(mn = mn0, std = std0, acptmn = acptmn, acptstd = acptstd)
  }
  
  drawpar_fast <- function(alpha0, beta0, theta0, acptalpha, acptbeta) {
    alpha1 <- alpha0 + cand_a * rnorm(J)
    beta1  <- beta0  + cand_b * rnorm(J)
    
    # current ll
    ll0 <- loglike_mat_fast(alpha0, beta0, theta0)
    
    # ---- alpha update (vectorized) ----
    ll_a1 <- loglike_mat_fast(alpha1, beta0, theta0)
    
    d0 <- colSums(ll0)   + logprior_alpha(alpha0)
    d1 <- colSums(ll_a1) + logprior_alpha(alpha1)
    
    acc <- d1 - d0
    accind <- (log(runif(J)) < acc)
    
    alpha0[accind] <- alpha1[accind]
    acptalpha[accind] <- acptalpha[accind] + 1L
    
    # update ll0 to reflect accepted alpha columns (avoid a full recompute)
    if (any(accind)) ll0[, accind] <- ll_a1[, accind]
    
    # ---- beta update (vectorized) ----
    ll_b1 <- loglike_mat_fast(alpha0, beta1, theta0)
    
    d0 <- colSums(ll0)   + logprior_beta(beta0)
    d1 <- colSums(ll_b1) + logprior_beta(beta1)
    
    acc <- d1 - d0
    accind <- (log(runif(J)) < acc)
    
    beta0[accind] <- beta1[accind]
    acptbeta[accind] <- acptbeta[accind] + 1L
    
    list(alpha = alpha0, beta = beta0, acpta = acptalpha, acptb = acptbeta)
  }
  
  # ---- State containers ----
  thstate  <- matrix(0, N, chains)
  astate   <- matrix(0, J, chains)
  bstate   <- matrix(0, J, chains)
  mnstate  <- matrix(0, ngrp, chains)
  stdstate <- matrix(0, ngrp, chains)
  
  sthchains   <- matrix(0, N, chains)
  sachains    <- matrix(0, J, chains)
  sbchains    <- matrix(0, J, chains)
  smnchains   <- matrix(0, ngrp, chains)
  sstdchains  <- matrix(0, ngrp, chains)
  
  sthchains2  <- matrix(0, N, chains)
  sachains2   <- matrix(0, J, chains)
  sbchains2   <- matrix(0, J, chains)
  smnchains2  <- matrix(0, ngrp, chains)
  sstdchains2 <- matrix(0, ngrp, chains)
  
  acptalpha <- rep(0L, J)
  acptbeta  <- rep(0L, J)
  acpttheta <- rep(0L, N)
  acptmn    <- 0L
  acptstd   <- 0L
  
  alphastates <- array(0, dim = c(J, (it - burn), chains))
  betastates  <- array(0, dim = c(J, (it - burn), chains))
  thetastates <- array(0, dim = c(N, (it - burn), chains))
  
  # ---- Initialize each chain ----
  for (c in 1:chains) {
    init <- initiate()
    astate[, c]   <- init$alpha
    bstate[, c]   <- init$beta
    thstate[, c]  <- init$theta
    mnstate[, c]  <- init$mn
    stdstate[, c] <- init$std
  }
  
  # ---- MCMC ----
  for (c in 1:chains) {
    for (iter in 1:it) {
      
      theta0 <- thstate[, c]
      alpha0 <- astate[, c]
      beta0  <- bstate[, c]
      mn0    <- mnstate[, c]
      std0   <- stdstate[, c]
      
      drawth <- drawtheta_fast(alpha0, beta0, theta0, mn0, std0, acpttheta)
      theta0 <- drawth$theta
      acpttheta <- drawth$acptth
      
      drawl <- drawlambda_fast(theta0, mn0, std0, acptmn, acptstd)
      mn0 <- drawl$mn
      std0 <- drawl$std
      acptmn <- drawl$acptmn
      acptstd <- drawl$acptstd
      
      drawp <- drawpar_fast(alpha0, beta0, theta0, acptalpha, acptbeta)
      alpha0 <- drawp$alpha
      beta0  <- drawp$beta
      acptalpha <- drawp$acpta
      acptbeta  <- drawp$acptb
      
      if (iter > burn) {
        idx <- iter - burn
        
        sthchains[, c]   <- sthchains[, c]   + theta0
        sthchains2[, c]  <- sthchains2[, c]  + theta0^2
        smnchains[, c]   <- smnchains[, c]   + mn0
        smnchains2[, c]  <- smnchains2[, c]  + mn0^2
        sstdchains[, c]  <- sstdchains[, c]  + std0
        sstdchains2[, c] <- sstdchains2[, c] + std0^2
        sachains[, c]    <- sachains[, c]    + alpha0
        sachains2[, c]   <- sachains2[, c]   + alpha0^2
        sbchains[, c]    <- sbchains[, c]    + beta0
        sbchains2[, c]   <- sbchains2[, c]   + beta0^2
        
        alphastates[, idx, c] <- alpha0
        betastates[,  idx, c] <- beta0
        thetastates[, idx, c] <- theta0
      }
      
      thstate[, c]  <- theta0
      astate[, c]   <- alpha0
      bstate[, c]   <- beta0
      mnstate[, c]  <- mn0
      stdstate[, c] <- std0
    }
  }
  
  # ---- Posterior means / SDs ----
  M <- (it - burn) * chains
  
  sthsum   <- rowSums(sthchains);   sthsum2   <- rowSums(sthchains2)
  sasum    <- rowSums(sachains);    sasum2    <- rowSums(sachains2)
  sbsum    <- rowSums(sbchains);    sbsum2    <- rowSums(sbchains2)
  smnsum   <- rowSums(smnchains);   smnsum2   <- rowSums(smnchains2)
  sstdsum  <- rowSums(sstdchains);  sstdsum2  <- rowSums(sstdchains2)
  
  alpha <- sasum / M
  beta  <- sbsum / M
  estpars <- cbind(alpha, beta)
  colnames(estpars) <- c("a-par", "b-par")
  
  a.sd <- sqrt((sasum2 - (sasum^2) / M) / ((it - burn - 1) * chains))
  b.sd <- sqrt((sbsum2 - (sbsum^2) / M) / ((it - burn - 1) * chains))
  sdpars <- cbind(a.sd, b.sd)
  colnames(sdpars) <- c("a-par", "b-par")
  
  theta <- sthsum / M
  th.sd <- sqrt((sthsum2 - (sthsum^2) / M) / ((it - burn - 1) * chains))
  
  mn <- smnsum / M
  mn.sd <- sqrt((smnsum2 - (smnsum^2) / M) / ((it - burn - 1) * chains))
  
  std <- sstdsum / M
  std.sd <- sqrt((sstdsum2 - (sstdsum^2) / M) / ((it - burn - 1) * chains))
  
  # ---- Acceptance ----
  AR.alpha <- acptalpha / (it * chains)
  AR.beta  <- acptbeta  / (it * chains)
  AR <- cbind(AR.alpha, AR.beta)
  colnames(AR) <- c("a-par", "b-par")
  
  AR.theta <- mean(acpttheta / (it * chains))
  AR.mn <- acptmn / (it * chains)
  AR.std <- acptstd / (it * chains)
  
  # ---- R-hat (same structure as your original) ----
  Rstat <- "NA; increase the number of chains to compute R stat"
  if (chains > 1) {
    Wa <- matrix(0, J, chains)
    Wb <- matrix(0, J, chains)
    for (c in 1:chains) {
      Sa <- (alphastates[, , c] - sachains[, c] / (it - burn))^2
      Sb <- (betastates[,  , c] - sbchains[, c] / (it - burn))^2
      Wa[, c] <- rowSums(Sa) / (it - burn - 1)
      Wb[, c] <- rowSums(Sb) / (it - burn - 1)
    }
    Walpha <- rowMeans(Wa)
    Wbeta  <- rowMeans(Wb)
    
    Ba <- matrix(0, J, chains)
    Bb <- matrix(0, J, chains)
    for (c in 1:chains) {
      Ba[, c] <- (sachains[, c] / (it - burn) - sasum / M)^2
      Bb[, c] <- (sbchains[, c] / (it - burn) - sbsum / M)^2
    }
    Balpha <- (it - burn) * rowSums(Ba) / (chains - 1)
    Bbeta  <- (it - burn) * rowSums(Bb) / (chains - 1)
    
    Var.a <- (1 - 1 / (it - burn)) * Walpha + (1 / (it - burn)) * Balpha
    Var.b <- (1 - 1 / (it - burn)) * Wbeta  + (1 / (it - burn)) * Bbeta
    
    Ralpha <- sqrt(Var.a / Walpha)
    Rbeta  <- sqrt(Var.b / Wbeta)
    
    Rstat <- cbind(Ralpha, Rbeta)
    colnames(Rstat) <- c("rhat a-par", "rhat b-par")
  }
  
  # ---- Model fit (log-likelihood at posterior mean parameters) ----
  ll <- sum(loglike_mat_fast(alpha, beta, theta))
  npar <- 2 * J
  AIC <- -2 * ll + npar
  BIC <- -2 * ll + npar * log(N)
  
  # DIC (kept as in spirit of your original; still expensive)
  Dbar <- 0
  for (c in 1:chains) {
    for (i in 1:(it - burn)) {
      Dbar <- Dbar + sum(loglike_mat_fast(alphastates[, i, c], betastates[, i, c], thetastates[, i, c]))
    }
  }
  Dbar <- -2 * Dbar / (chains * (it - burn))
  pD <- Dbar - (-2 * ll)
  DIC <- Dbar + pD
  model.fit <- cbind(AIC, BIC, DIC)
  
  return(list(
    item.par    = estpars,
    item.sd     = sdpars,
    item.accept = AR,
    item.conv   = Rstat,
    theta.est   = theta,
    theta.psd   = th.sd,
    theta.accept= AR.theta,
    mn.est      = mn,
    mn.psd      = mn.sd,
    mn.accept   = AR.mn,
    SD.est      = std,
    SD.psd      = std.sd,
    SD.accept   = AR.std,
    model.fit   = model.fit,
    a.sample    = alphastates,
    b.sample    = betastates,
    th.sample   = thetastates
  ))
}
