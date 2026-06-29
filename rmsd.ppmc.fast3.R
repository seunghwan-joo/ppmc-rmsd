rmsd.ppmc <- function(data, mn = 0, std = 1,
                      a.sample,
                      b.sample,
                      th.sample) {
  
  # ----------------------------
  # basic setup
  # ----------------------------
  X <- as.matrix(data)
  itemnames <- colnames(X)
  N <- nrow(X)
  J <- ncol(X)
  
  alphastates <- as.matrix(a.sample)
  betastates  <- as.matrix(b.sample)
  thetastates <- as.matrix(th.sample)
  
  n.iter <- ncol(alphastates)
  
  if (nrow(alphastates) != J) stop("a.sample must have J rows.")
  if (nrow(betastates)  != J) stop("b.sample must have J rows.")
  if (nrow(thetastates) != N) stop("th.sample must have N rows.")
  if (ncol(betastates)  != n.iter) stop("a.sample and b.sample must have same number of columns.")
  if (ncol(thetastates) != n.iter) stop("a.sample and th.sample must have same number of columns.")
  
  # ----------------------------
  # 2PL probability function
  # ----------------------------
  twopl <- function(a, b, th) {
    eta <- outer(th, a, "*") - matrix(a * b,
                                      nrow = length(th),
                                      ncol = length(a),
                                      byrow = TRUE)
    plogis(eta)
  }
  
  # ----------------------------
  # quadrature setup
  # ----------------------------
  qt <- seq(-4, 4, by = 0.2)
  n.qt <- length(qt)
  
  grpwt <- dnorm(qt, mean = mn, sd = std)
  grpwt <- grpwt / sum(grpwt)
  log_grpwt <- log(grpwt)
  
  range.ind <- which(round(grpwt, 2) %in% 0.01)
  minind <- range.ind[1]
  maxind <- range.ind[length(range.ind)]
  use_idx <- minind:maxind
  
  # ----------------------------
  # precompute observed-data masks
  # ----------------------------
  obs_mask <- !is.na(X)              # N x J logical
  X_num <- X
  X_num[!obs_mask] <- 0              # replace NA by 0 for matrix ops
  
  X1_obs <- (X_num == 1) * 1         # indicator for observed 1
  X0_obs <- (X_num == 0 & obs_mask) * 1   # indicator for observed 0
  
  # ----------------------------
  # helper: compute RMSD from a response matrix
  # ----------------------------
  compute_rmsd <- function(X1, X0, X_num_use, obs_mask_use, exp.icc, logp, logq) {
    
    # log-likelihood for all persons over all quadrature points
    # N x n.qt
    loglik <- X1 %*% t(logp) + X0 %*% t(logq)
    
    # posterior in log-space
    logpost <- sweep(loglik, 2, log_grpwt, "+")
    
    # stabilize row-wise
    row_max <- apply(logpost, 1, max)
    post <- exp(logpost - row_max)
    post <- post / rowSums(post)
    
    # n.qt x J
    pseudo <- t(post) %*% X_num_use
    deno   <- t(post) %*% obs_mask_use
    
    obs.icc <- pseudo / deno
    
    rmsd <- sqrt(colSums((obs.icc[use_idx, , drop = FALSE] -
                            exp.icc[use_idx, , drop = FALSE])^2 *
                           grpwt[use_idx]))
    rmsd
  }
  
  # ----------------------------
  # storage
  # ----------------------------
  sample.obs.rmsd <- matrix(NA_real_, n.iter, J)
  sample.rep.rmsd <- matrix(NA_real_, n.iter, J)
  
  # ----------------------------
  # main loop over posterior draws
  # ----------------------------
  for (rep in seq_len(n.iter)) {
    
    rep.alpha <- alphastates[, rep]
    rep.beta  <- betastates[, rep]
    rep.th    <- thetastates[, rep]
    
    # probabilities for replicated data generation
    rep.PRval <- twopl(rep.alpha, rep.beta, rep.th)   # N x J
    
    # generate replicated data (vectorized)
    rep.X <- (matrix(runif(N * J), N, J) < rep.PRval) * 1
    
    # ICC at quadrature points
    exp.icc <- twopl(rep.alpha, rep.beta, qt)         # n.qt x J
    
    # log probabilities for observed/replicated posterior computations
    # clamp to avoid log(0)
    eps <- 1e-12
    exp.icc2 <- pmin(pmax(exp.icc, eps), 1 - eps)
    logp <- log(exp.icc2)
    logq <- log1p(-exp.icc2)
    
    # observed RMSD
    obs.rmsd <- compute_rmsd(
      X1 = X1_obs,
      X0 = X0_obs,
      X_num_use = X_num,
      obs_mask_use = obs_mask * 1,
      exp.icc = exp.icc,
      logp = logp,
      logq = logq
    )
    
    # replicated RMSD
    X1_rep <- (rep.X == 1) * 1
    X0_rep <- (rep.X == 0) * 1
    
    rep.rmsd <- compute_rmsd(
      X1 = X1_rep,
      X0 = X0_rep,
      X_num_use = rep.X,
      obs_mask_use = matrix(1, N, J),
      exp.icc = exp.icc,
      logp = logp,
      logq = logq
    )
    
    sample.obs.rmsd[rep, ] <- obs.rmsd
    sample.rep.rmsd[rep, ] <- rep.rmsd
  }
  
  # posterior predictive p-values
  ppp.rmsd <- matrix(round(colMeans(sample.rep.rmsd > sample.obs.rmsd, na.rm = TRUE), 2),
                     ncol = 1)
  rownames(ppp.rmsd) <- itemnames
  colnames(ppp.rmsd) <- "RMSD PPP-value"
  
  colnames(sample.obs.rmsd) <- itemnames
  colnames(sample.rep.rmsd) <- itemnames
  
  return(list(
    obs.rmsd = sample.obs.rmsd,
    rep.rmsd = sample.rep.rmsd,
    ppp = ppp.rmsd
  ))
}